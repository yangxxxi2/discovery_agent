from storage.db import DB
from dotenv import load_dotenv
from knowledge_base_constructor import KnowledgeBaseConstructor
from question_analogical_augmenter import ResearchQuestionAnalogicalAugmenter
from response_generator import ResponseGenerator

load_dotenv()

if __name__ == "__main__":
    
    # ============================================================
    # 阶段1: 构建知识库（从已发表的研究中学习）
    # ============================================================

    db = DB()
    db.reset_db()

    knowledge_base = KnowledgeBaseConstructor()
    abstracts = knowledge_base.fetch_abstracts_from_pubmed(max_results=10)

    success_count = 0
    for item in abstracts:
        if(knowledge_base.persist_parse_abstract(item, db)):
            success_count += 1
    print(f"\n✓ 成功处理 {success_count} 篇摘要，知识库构建完成\n")

    # ============================================================
    # 阶段2: 用户问题增强（举一反三生成创新问题）
    # ============================================================
    augmenter = ResearchQuestionAnalogicalAugmenter()
    generator = ResponseGenerator()
    
    user_query = "What is the effectiveness of cognitive behavioral therapy in treating adolescent depression?"
    
    print("→ 步骤1: 解析用户问题...")
    structured_user_query = augmenter.parse_user_query(user_query, db)
    
    print("\n→ 步骤2: 检查问题模型完整性...")
    is_completed = augmenter.check_pending_model_slot_filled(structured_user_query)
    
    if is_completed:
        print("\n✓ 最终: 问题模型已完整，直接生成评估\n")
        generator.generate_question_evaluation(structured_user_query)
    else:
        print("\n→ 步骤3: 生成平行模型（多角度查询）...")
        variant_models = augmenter.generate_variant_models_by_extend_filled_slot(structured_user_query)
        print(f"  生成了 {len(variant_models)} 个平行查询模型")
        
        print("\n→ 步骤4: RAG检索 - 学习已有研究的结构化模式...")
        evidence = augmenter.retrieve_related_evidence_based_variant_models(variant_models, db, 10, 0.5)
        print(f"\n  检索到 {len(evidence)} 个相关研究作为学习样本, 相似度最高的3个:")
        for i, e in enumerate(evidence[:3], 1):
            print(f"    {i}. {e['title'][:60]}... (相似度: {e['similarity_score']:.2f})")
            if 'components' in e:
                print(f"       结构: {e['components']}")
        
        print("\n→ 步骤5: 举一反三生成创新性研究问题...")
        print("  基于学习到的模式:")
        print("    1) 对比分析: 现有研究用了什么组合?")
        print("    2) Gap识别: 哪些方向未被充分研究?")
        print("    3) 创新生成: 提出有价值的新研究问题")
        candidate_models = augmenter.augment_candidate_models_by_analogical_generation(structured_user_query, evidence)
        print(f"  生成了 {len(candidate_models)} 个候选研究问题")
        
        print("\n→ 步骤6: 评估和排序候选问题...")
        ranked_candidate_models = augmenter.evaluate_and_rank_candidate_models(candidate_models)
        print(f"  完成评估，排序如下:")
        for i, ranked in enumerate(ranked_candidate_models, 1):
            model = ranked['model']
            score = ranked.get('total_score', 0)
            print(f"\n  排名 {i} (总分: {score}):")
            if 'innovation_rationale' in model:
                print(f"    创新点: {model['innovation_rationale']}")
        
        print("\n✓ 最终: 生成最终研究问题...")
        generator.generate_research_questions(ranked_candidate_models)
        
    db.close()