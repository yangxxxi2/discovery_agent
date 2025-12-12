import sys
sys.path.append('/Users/xinyang/PythonProjects/discovery_agent')

from question_analogical_augmenter import ResearchQuestionAnalogicalAugmenter
from storage.db import DB

def test_generate_parallel_models():
    """测试生成平行模型功能"""
    print("\n=== 测试 generate_parallel_models_by_extend_filled_slot ===\n")
    
    db = DB()
    augmenter = ResearchQuestionAnalogicalAugmenter()
    
    # 模拟一个部分填充的 PICO 模型
    pending_model = {
        "user_query": "青少年抑郁症认知行为疗法的效果",
        "task_type": "generate_question",
        "research_type": "quantitative",
        "model_type": "PICO",
        "P": "青少年抑郁症患者",
        "I": "认知行为疗法",
        "C": "Unknown",
        "O": "Unknown"
    }
    
    print("输入的待补全模型:")
    print(f"  P: {pending_model['P']}")
    print(f"  I: {pending_model['I']}")
    print(f"  C: {pending_model['C']}")
    print(f"  O: {pending_model['O']}")
    print()
    
    # 生成平行模型
    parallel_models = augmenter.generate_variant_models_by_extend_filled_slot(
        pending_model=pending_model,
        db=db
    )
    
    print(f"\n生成的平行模型数量: {len(parallel_models)}")
    print("\n平行模型详情:")
    for i, model in enumerate(parallel_models, 1):
        print(f"\n模型 {i}:")
        for key, value in model.items():
            if key != "model_type":
                print(f"  {key}: {value}")
    
    db.close()
    return parallel_models


def test_full_workflow():
    """测试完整的增强工作流"""
    print("\n=== 测试完整的研究问题增强工作流 ===\n")
    
    db = DB()
    augmenter = ResearchQuestionAnalogicalAugmenter()
    
    # 用户查询
    user_query = "青少年抑郁症认知行为疗法的疗效如何？"
    print(f"用户输入: {user_query}\n")
    
    # 步骤1: 解析用户查询
    print("步骤1: 解析用户查询...")
    pending_model = augmenter.parse_user_query(user_query, db)
    print(f"  任务类型: {pending_model.get('task_type')}")
    print(f"  研究类型: {pending_model.get('research_type')}")
    print(f"  模型类型: {pending_model.get('model_type')}")
    print()
    
    # 步骤2: 检查槽位填充情况
    print("步骤2: 检查槽位填充...")
    is_complete = augmenter.check_pending_model_slot_filled(pending_model)
    print(f"  模型完整: {is_complete}\n")
    
    if not is_complete:
        # 步骤3: 生成平行模型
        print("步骤3: 生成平行模型（扩展已填充槽位）...")
        parallel_models = augmenter.generate_variant_models_by_extend_filled_slot(
            pending_model=pending_model,
            db=db
        )
        print()
        
        # 步骤4: RAG检索相关证据
        print("步骤4: RAG检索相关证据...")
        try:
            evidence = augmenter.retrieve_related_evidence_based_variant_models(
                variant_models=parallel_models,
                db=db,
                top_k=5,
                distance_threshold=0.8
            )
            print(f"\n检索到的证据数量: {len(evidence)}")
            
            if evidence:
                print("\n前3条证据:")
                for i, ev in enumerate(evidence[:3], 1):
                    print(f"\n证据 {i}:")
                    print(f"  相似度: {ev['similarity_score']:.3f}")
                    print(f"  研究类型: {ev['research_type']}")
                    print(f"  标题: {ev['title'][:80]}...")
                    print(f"  组件: {ev['components']}")
            else:
                print("  注意: 数据库中暂无相关研究，无法演示RAG检索")
                print("  建议先使用 knowledge_base_constructor.py 构建知识库")
        except Exception as e:
            print(f"  RAG检索出错: {e}")
            print("  可能原因: pgvector 未安装或数据库中没有向量数据")
    
    db.close()
    print("\n测试完成！")


if __name__ == "__main__":
    # 测试1: 平行模型生成
    test_generate_parallel_models()
    
    print("\n" + "="*60 + "\n")
    
    # 测试2: 完整工作流
    test_full_workflow()
