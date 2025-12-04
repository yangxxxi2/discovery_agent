import sys
import os
# 添加项目根目录到路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from storage import DB
from dotenv import load_dotenv

# 加载环境变量
load_dotenv()

db = DB()

class TestDB:

    def test_reset_db():
        print(f"数据库配置: {db.config}")
        db.reset_db()
    
    def test_get_frameworks_and_elements():
        frameworks = db.get_all_frameworks()
        print("已加载的框架:")
        for fw in frameworks:
            print(f"- {fw['label']}: {fw['description']}")
            elements = db.get_framework_elements(fw['id'])
            for elem in elements:
                print(f"  [{elem['label']}] {elem['description']}")

    def test_insert_question_source():
        db.insert_question_source(
            source_label='PubMed',
            source_id='PMID12345678',
            research_type='Randomized Controlled Trials',
            grade='A',
            title='Effect of Intervention X on Outcome Y in Population Z',
            abstract='This randomized controlled trial investigated the effect of intervention X...'
        )
    
    def test_get_question_source():
        db.insert_question_source(
            source_label='ClinicalTrials.gov',
            source_id='ased1',
            research_type='RWS',
            grade='A',
            title='aaaa',
            abstract='xxxxx'
        )
        source = db.get_question_source_by_title('ClinicalTrials.gov', 'aaaa')
        print(source)

    def test_insert_question():
        import numpy as np
        question_source = db.get_question_source_by_title('ClinicalTrials.gov', 'aaaa')
        framework1 = db.get_framework_by_label('PICO')
        framework2 = db.get_framework_by_label('SPIDER')
        if question_source and framework1:
            embedding = np.random.rand(1536).tolist()
            db.insert_question(
                question_source_id=question_source['id'],
                framework_id=framework1['id'],
                embedding=embedding
            )
        if question_source and framework2:
            db.insert_question(
                question_source_id=question_source['id'],
                framework_id=framework2['id']
            )

    def test_insert_question_element():
        pico_framework = db.get_framework_by_label('PICO')
        pico_elements = db.get_framework_elements(pico_framework['id'])
        pico_values = {
            'P': 'Adults with type 2 diabetes aged 40-65',
            'I': 'Metformin 1000mg twice daily',
            'C': 'Placebo',
            'O': 'HbA1c reduction after 12 weeks'
        }
        
        for element in pico_elements:
            if element['label'] in pico_values:
                db.insert_question_element(
                    question_id=1,
                    framework_element_id=element['id'],
                    question_element_label=pico_values[element['label']]
                )

    def test_get_question_elements_by_question_id():
        question_elements = db.get_question_elements_by_question_id(1)
        print("Question Elements for Question ID 1:")
        for elem in question_elements:
            print(f"- [{elem['question_element_label']}]")

    def test_get_statistics():
        stats = db.get_statistics()
        print("数据库统计信息:")
        for key, value in stats.items():
            print(f"{key}: {value}")

if __name__ == "__main__":
    test_db = TestDB()
    test_db.test_reset_db()
    # test_db.test_get_frameworks_and_elements()
    # test_db.test_insert_question_source()
    # test_db.test_get_question_source()
    # test_db.test_insert_question()
    # test_db.test_insert_question_element()
    # test_db.test_get_question_elements_by_question_id()
    # test_db.test_get_statistics()