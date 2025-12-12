from typing import List
from .question_model_extractor import QuestionModelExtractor
from storage.db import DB
from .llm import llm
import json
from .prompt import EXTEND_FILLED_SLOT_PROMPT, ANALOGICAL_GENERATION_PROMPT, EVALUATE_CANDIDATE_MODELS_PROMPT

class ResearchQuestionAnalogicalAugmenter:
    """
    研究问题增强器 - 通过举一反三生成创新性研究问题
    
    核心思想：
    1. 从知识库检索相关研究，学习其模式
    2. 识别研究gap和未探索的方向
    3. 生成创新性的、有价值的研究问题（而非重复现有研究）
    """

    def __init__(self):
        self.extractor = QuestionModelExtractor()

    def parse_user_query(self, user_query: str, db: DB) -> dict:
        task_type = self.extractor.determine_task_type(user_query)
        question_model = self.extractor.extract_question_model(user_query, task_type, db) if task_type != None else {}
        research_type = self.extractor.extract_research_type(user_query) if task_type != None else ""
        return {
            "user_query": user_query,
            "task_type": task_type,
            "research_type": research_type,
            **question_model,
        }

    def check_pending_model_slot_filled(self, pending_model: dict) -> bool:
        for key, value in pending_model.items():
            if key in ["user_query", "task_type", "research_type", "model_type"]:
                continue
            if value == "Unknown" or value == "" or value is None:
                print(f"未填充的位置: {key}")
                return False
        print("所有位置均已填充")
        return True
    
    def generate_variant_models_by_extend_filled_slots(self, pending_model: dict) -> List[dict]:
        """
        通过扩展已填充槽位生成平行模型
        
        策略：
        1. 识别所有已填充的槽位
        2. 使用LLM为每个槽位生成2-3个语义相似的变体
        3. 组合生成多个平行的研究模型
        
        示例：
        原始: P="青少年抑郁症"
        变体: ["青春期情绪障碍", "青年抑郁症状", "未成年人抑郁"]
        
        目的：扩大RAG检索范围，发现更多相关研究
        """
        parallel_models = []
        base_model = {k: v for k, v in pending_model.items() 
                     if k not in ["user_query", "task_type", "research_type"]}
        parallel_models.append(base_model)

        filled_slots = {}
        for key, value in pending_model.items():
            if key not in ["user_query", "task_type", "research_type", "model_type"]:
                if value and value != "Unknown" and value != "":
                    filled_slots[key] = value
        
        if not filled_slots:
            return parallel_models
        
        formatted_prompt = EXTEND_FILLED_SLOTS_PROMPT.format()
       
        try:
            result = llm.invoke(formatted_prompt)
            result = result.content if hasattr(result, 'content') else json.loads(str(result))
     
            for slot_name, variant_list in result.items():
                if slot_name in filled_slots:
                    for variant_value in variant_list:
                        variant_model = base_model.copy()
                        variant_model[slot_name] = variant_value
                        parallel_models.append(variant_model)
            
            print(f"✓ 生成了 {len(parallel_models)} 个平行模型（包含原始模型）")
            return parallel_models
            
        except Exception as e:
            print(f"生成槽位变体时出错: {e}")
            return parallel_models

    def retrieve_related_evidence_by_variant_models(self, variant_models: List[dict], db: DB, top_k: int, distance_threshold: float) -> List[dict]:
        """
        基于平行模型检索相关证据 (RAG检索)
        
        Args:
            models: 平行模型列表
            db: 数据库连接
            top_k: 每个模型检索的最大结果数
            distance_threshold: 距离阈值（< 0.8 表示高相关）
        
        Returns:
            相关证据列表，按相似度排序
        """
        evidence_list = []
        
        for variant in variant_models:
            try:
                model_embedding = self.extractor.embed_question_model(variant)
                if len(model_embedding) == 0:
                    continue
                
                similar_questions = db.search_similar_questions(
                    model_embedding=model_embedding,
                    limit=top_k,
                    distance_threshold=distance_threshold
                )
                
                for q in similar_questions:
                    question_elements = db.get_question_elements_by_question_id(q['id'])
                    question_source = db.get_question_source_by_id(q['question_source_id'])
                    framework = db.get_framework_by_id(q['framework_id'])
                    
                    evidence = {
                        'question_id': q['id'],
                        'model_type': framework['label'] if framework else 'Unknown',
                        'similarity_score': 1 - q['distance'],  # 转换为相似度
                        'research_type': question_source.get('research_type', '') if question_source else '',
                        'title': question_source.get('title', '') if question_source else '',
                        'abstract': question_source.get('abstract', '') if question_source else '',
                        'components': {}
                    }
                    
                    for elem in question_elements:
                        elem_label = elem.get('element_label') or elem.get('label', 'Unknown')
                        evidence['components'][elem_label] = elem.get('question_element_label', '')
                    
                    evidence_list.append(evidence)
                    
            except Exception as e:
                print(f"处理模型时出错: {e}")
                continue

        seen_ids = set()
        unique_evidence = []
        for ev in sorted(evidence_list, key=lambda x: x['similarity_score'], reverse=True):
            if ev['question_id'] not in seen_ids:
                seen_ids.add(ev['question_id'])
                unique_evidence.append(ev)
        
        print(f"✓ RAG检索完成：从 {len(variant_models)} 个平行模型中检索到 {len(unique_evidence)} 条独特证据")
        return unique_evidence[:top_k]

    def fill_pending_model_by_compose_evidence_components(self, evidence: List[dict]) -> dict:
       #TODO:识别unfilled slots；从evidence提取相应slot对应的component，并构建pool；智能采样组合，填充pending model
        pass

    def augment_candidate_models_by_analogical_generation(self, pending_model: dict, evidence: List[dict]) -> List[dict]:
        #TODO：gap模板驱动的举一反三式生成
        #构建受控算子库；对每个composed model识别gap，计算gap score；重复的研究则进行受控调整，新颖的研究则保留
        """
        生成完整模型 - 通过举一反三补全未填充槽位
        
        核心创新点：相比传统RAG，不是直接照搬现有研究，而是分析已有研究的pattern和gap，生成创新性的、有研究价值的研究问题
        
        策略：
        - Pattern分析：相似研究用了什么P/I/C/O组合？
        - Gap识别：哪些组合还没有被研究？
        - 创新生成：提出有价值的新组合
        """


        fromatted_prompt = EXTEND_FILLED_SLOTS_PROMPT.format()
        try:
            result = llm.invoke(fromatted_prompt)
            result = result.content if hasattr(result, 'content') else json.load(str(result))

            for model in result:
                model['user_query'] = pending_model['user_query']
                model['task_type'] = pending_model['task_type']
                model['research_type'] = pending_model.get('research_type', 'Unknown')
            
            return result
            
        except Exception as e:
            print(f"生成候选模型时出错: {e}")
            return []

    def evaluate_and_rank_candidate_models(self, candidate_models: List[dict]) -> List[dict]:
        """
        评估和排序生成的候选模型
        
        评估维度：
        1. Novelty：与现有研究的差异度
        2. Feasibility：研究的可执行性
        3. Interesting：研究问题的吸引力
        4. Relevant：与用户查询的相关性
        5. Ethical Risk：是否符合研究伦理
        """
        
        fromatted_prompt = EVALUATE_CANDIDATE_MODELS_PROMPT.format()
        try:
            result = llm.invoke(fromatted_prompt)
            result = result.content if hasattr(result, 'content') else json.load(str(result))
          
            return result
            
        except Exception as e:
            print(f"评估候选模型时出错: {e}")
            return [{'model': m, 'total_score': 0, 'rationale': 'Evaluation failed'} 
                    for m in candidate_models]
