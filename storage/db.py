import os
import psycopg2
from psycopg2.extras import RealDictCursor
import logging
from typing import List, Dict, Optional
from ..constants import (
    QUESTION_TYPE_INTEVENTION_TIME,
    QUESTION_TYPE_INTEVENTION,
    QUESTION_TYPE_OBSERVATION_TIME,
    QUESTION_TYPE_OBSERVATION,
    QUESTION_TYPE_PREVALENCE,
    QUESTION_TYPE_DIAGNOSIS,
    QUESTION_TYPE_QUALITATION,
    QUESTION_TYPE_HEALTH_PRACTICE,
    QUESTION_TYPE_HEALTH_SERVICE,
    QUESTION_MODEL_PICO,
    QUESTION_MODEL_PICOT,
    QUESTION_MODEL_PECO,
    QUESTION_MODEL_PEO,
    QUESTION_MODEL_CoCoPop,
    QUESTION_MODEL_PIRD,
    QUESTION_MODEL_SPIDER,
    QUESTION_MODEL_SPICE,
    QUESTION_MODEL_ECLIPSe,
)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DB:

    def __init__(self):
        self.config = {
            "host": os.getenv('DB_HOST'),
            "port": int(os.getenv('DB_PORT')),
            "database": os.getenv('DB_NAME'),
            "user": os.getenv('DB_USER'),
            "password": os.getenv('DB_PASSWORD')
        }
        self.connection = psycopg2.connect(**self.config)
    
    def reset_db(self):
        self.drop_tables()
        self.create_tables()
        self.initialize_default_frameworks()
    
    def commit(self):
        try:
            self.connection.commit()
        except Exception as e:
            if self.connection:
                self.connection.rollback()
            logger.error(f"数据库操作失败: {e}")
            raise
    
    def close(self):
        if self.connection:
            self.connection.close()

    def create_tables(self):
        create_query = """
        CREATE EXTENSION IF NOT EXISTS vector;

        CREATE TABLE IF NOT EXISTS framework (
            id                      SERIAL PRIMARY KEY,
            label                   VARCHAR(50) UNIQUE NOT NULL,
            question_type           VARCHAR(100) NOT NULL,
            description             TEXT,
            example                 TEXT,
            CONSTRAINT unique_framework UNIQUE(label)
        );

        CREATE TABLE IF NOT EXISTS framework_elements (
            id                      SERIAL PRIMARY KEY,
            framework_id            INT NOT NULL REFERENCES framework(id) ON DELETE CASCADE,
            label                   VARCHAR(50) NOT NULL,
            "order"                 INT NOT NULL,
            description             TEXT,
            CONSTRAINT unique_framework_element UNIQUE(framework_id, label)
        );

        CREATE TABLE IF NOT EXISTS question_source(
            id                      SERIAL PRIMARY KEY,
            source_label            VARCHAR(50) NOT NULL,
            source_id               VARCHAR(100),
            research_type           VARCHAR(100) NOT NULL,
            question_type           VARCHAR(100) NOT NULL,
            grade_level             VARCHAR(50),
            ocebm_2011_level        VARCHAR(50),
            title                   TEXT NOT NULL,
            abstract                TEXT NOT NULL,
            CONSTRAINT unique_source UNIQUE(source_label, title)
        );

        CREATE TABLE IF NOT EXISTS question (
            id                      SERIAL PRIMARY KEY,
            question_source_id      INT NOT NULL REFERENCES question_source(id) ON DELETE CASCADE,
            framework_id            INT NOT NULL REFERENCES framework(id),
            embedding               VECTOR(1536),
            CONSTRAINT unique_question UNIQUE(question_source_id, framework_id)
        );

        CREATE TABLE IF NOT EXISTS question_elements (
            id                      SERIAL PRIMARY KEY,
            question_id             INT NOT NULL REFERENCES question(id) ON DELETE CASCADE,
            framework_element_id    INT NOT NULL REFERENCES framework_elements(id),
            question_element_label  TEXT NOT NULL,
            CONSTRAINT unique_question_element UNIQUE(question_id, framework_element_id)
        );

        CREATE INDEX IF NOT EXISTS idx_question_embedding_vector 
        ON question USING ivfflat (embedding vector_cosine_ops) 
        WITH (lists = 100);
        """
        with self.connection.cursor() as cur:
            cur.execute(create_query)
        self.commit()

    def drop_tables(self):
        tables = [
            'question_elements',
            'question',
            'question_source',
            'framework_elements',
            'framework'
        ]
        with self.connection.cursor() as cur:
            for table in tables:
                cur.execute(f"DROP TABLE IF EXISTS {table} CASCADE")
        self.commit()
        logger.info("所有数据表已删除")

    def initialize_default_frameworks(self) -> None:
        insert_framework_query = """
                        INSERT INTO framework (label, question_type, description, example)
                        VALUES (%s, %s, %s, %s)
                        ON CONFLICT (label) DO NOTHING
                        RETURNING id
                        """
        insert_framework_element_query = """
                        INSERT INTO framework_elements (framework_id, label, "order", description)
                        VALUES (%s, %s, %s, %s)
                        ON CONFLICT (framework_id, label) DO NOTHING
                        """
        search_framework_query = "SELECT id FROM framework WHERE label = %s"
        frameworks = [
            {
                'label': QUESTION_MODEL_PICOT,
                'question_type': QUESTION_TYPE_INTEVENTION_TIME,
                'description': '干预与结局关联分析，有时间维度',
                'example': '(1)在绝经后骨质疏松女性群体中，接受地诺单抗治疗相比安慰剂，在24个月内，是否能降低椎体骨折发生率;',
                'elements': [
                    {'label': 'P', 'order': 1, 'description': 'Population - 研究人群'},
                    {'label': 'I', 'order': 2, 'description': 'Intervention - 干预措施'},
                    {'label': 'C', 'order': 3, 'description': 'Comparison - 对照组'},
                    {'label': 'O', 'order': 4, 'description': 'Outcome - 结局指标'},
                    {'label': 'T', 'order': 5, 'description': 'Timeframe - 时间'}
                ]
            },
            {
                'label': QUESTION_MODEL_PICO,
                'question_type': QUESTION_TYPE_INTEVENTION,
                'description': '干预与结局关联分析，没有时间维度',
                'example': '(1)在2型糖尿病成人群体中，相比单独运动管理，GLP-1受体激动剂治疗是否能改善HbA1c并减少体重;',
                'elements': [
                    {'label': 'P', 'order': 1, 'description': 'Population - 研究人群'},
                    {'label': 'I', 'order': 2, 'description': 'Intervention - 干预措施'},
                    {'label': 'C', 'order': 3, 'description': 'Comparison - 对照组'},
                    {'label': 'O', 'order': 4, 'description': 'Outcome - 结局指标'}
                ]
            },
            {
                'label': QUESTION_MODEL_PEO,
                'question_type': QUESTION_TYPE_OBSERVATION,
                'description': '暴露与结局关联分析，没有对照组',
                'example': '(1)在暴露于二手烟的儿童群体中，哮喘发生率是否增加;(2)孕妇孕期生活在空气污染环境中与早产发生是否相关;',
                'elements': [
                    {'label': 'P', 'order': 1, 'description': 'Population - 研究人群'},
                    {'label': 'E', 'order': 2, 'description': 'Exposure - 暴露因素'},
                    {'label': 'O', 'order': 3, 'description': 'Outcome - 结局指标'}
                ]
            },
            {
                'label': QUESTION_MODEL_PECO,
                'question_type': QUESTION_TYPE_OBSERVATION_TIME,
                'description': '暴露与结局关联分析，有对照组',
                'example': '(1) 在吸烟成人群体中，相比不吸烟者，冠心病发生率是否增加;(2)孕妇孕期暴露于高浓度空气污染环境中，与低出生体重发生是否相关;',
                'elements': [
                    {'label': 'P', 'order': 1, 'description': 'Population - 研究人群'},
                    {'label': 'E', 'order': 2, 'description': 'Exposure - 暴露因素'},
                    {'label': 'C', 'order': 3, 'description': 'Comparison - 对照组'},
                    {'label': 'O', 'order': 4, 'description': 'Outcome - 结局指标'}
                ]
            },
            {
                'label': QUESTION_MODEL_CoCoPop,
                'question_type': QUESTION_TYPE_PREVALENCE,
                'description': '发病率/患病率的描述性分析；人群层面暴露与结局关联分析',
                'example': '(1)在40-70岁社区居民中，城市社区基层卫生服务中心体检环境下的高血压的患病率;'+
                           '(2)高中学生群体中抑郁症状在学校环境下的发生率;'+
                           '(3)在ICU住院患者中，医院获得性感染在重症监护环境下的发生率/患病率是多少',
                'elements': [
                    {'label': 'Co', 'order': 1, 'description': 'Condition - 疾病/健康问题'},
                    {'label': 'Co', 'order': 2, 'description': 'Context - 研究背景/场所/地区'},
                    {'label': 'P', 'order': 3, 'description': 'Population - 研究人群'}
                ]
            },
            {
                'label': QUESTION_MODEL_PIRD,
                'question_type': QUESTION_TYPE_DIAGNOSIS,
                'description': '诊断方法准确性评价',
                'example': '(1)高敏肌钙蛋白T（hs-cTnT）的0/1小时算法与冠脉造影方法相比，在诊断急诊胸痛成人患者的急性心肌梗死的表现;'+
                           '(2)超声弹性成像技术与肝活检方法相比，在诊断慢性肝病患者肝纤维化程度的准确性;'+
                           '(3)在疑似乳腺癌的女性中，乳腺MRI相比术后病理诊断，在诊断乳腺恶性肿瘤的准确性;',
                'elements': [
                    {'label': 'P', 'order': 1, 'description': 'Population - 研究人群'},
                    {'label': 'I', 'order': 2, 'description': 'Index test - 诊断方法'},
                    {'label': 'R', 'order': 3, 'description': 'Reference test - 金标准方法'},
                    {'label': 'D', 'order': 4, 'description': 'Diagnosis of interest - 诊断目标'}
                ]
            },
            {
                'label': QUESTION_MODEL_SPIDER,
                'question_type': QUESTION_TYPE_QUALITATION,
                'description': '关注患者体验感受、医生观点、医生经验等访谈/观察类质性研究，其目的在于探索现象本质，理解主观体验',
                'example': '(1)开展一项定性研究（r），对接受髋关节置换术的老年患者（s）进行半结构式访谈（d），主题分析患者对术后疼痛管理（pi）的体验感受（e）;' +
                           '(2)探索2型糖尿病患者是如何体验与理解日常血糖自我管理行为，采用半结构式访谈收集他们对自我管理挑战与促进因素的看法，深入了解其主观体验;',
                'elements': [
                    {'label': 'S', 'order': 1, 'description': 'Sample - 样本/研究对象'},
                    {'label': 'PI', 'order': 2, 'description': 'Phenomenon of Interest - 关注现象'},
                    {'label': 'D', 'order': 3, 'description': 'Design - 研究设计'},
                    {'label': 'E', 'order': 4, 'description': 'Evaluation - 评估'},
                    {'label': 'R', 'order': 5, 'description': 'Research type - 研究类型'}
                ]
            },
            {
                'label': QUESTION_MODEL_SPICE,
                'question_type': QUESTION_TYPE_HEALTH_PRACTICE,
                'description': '聚焦医疗健康领域实践服务的应用型研究，其目的在于评估、优化健康实践服务，提升医疗质量、患者结局与资源利用率，例如教育、护理流程、随访模式的实践评价（服务质量、满意度、依从性、流程效率等）',
                'example': '(1)在急诊科环境中（s），对于成人急诊患者（p），实施快速评估与分诊流程（i）相比传统护理流程（c），能否提高患者满意度和缩短等待时间（e）;'+
                           '(2)在社区卫生服务中心环境中，对于2型糖尿病患者，实施个体化健康教育干预相比常规健康教育，能否改善患者的血糖控制和生活质量;'+
                           '(3)在三级医院风湿免疫科，对于系统性红斑狼疮患者，实施护士主导的疾病管理门诊，相比传统医生门诊随访，能否改善疾病活动度控制与复诊依从性;',
                'elements': [
                    {'label': 'S', 'order': 1, 'description': 'Setting - 环境'},
                    {'label': 'P', 'order': 2, 'description': 'Perspective - 对象'},
                    {'label': 'I', 'order': 3, 'description': 'Intervention - 干预'},
                    {'label': 'C', 'order': 4, 'description': 'Comparison - 对照'},
                    {'label': 'E', 'order': 5, 'description': 'Evaluation - 评估'}
                ]
            },
            {
                'label': QUESTION_MODEL_ECLIPSe,
                'question_type': QUESTION_TYPE_HEALTH_SERVICE,
                'description': '聚焦卫生服务体系、功能、政策的交叉学科研究，其目的在于分析卫生服务的供给、需求、利用与分配，优化资源配置、完善卫生政策、提高公平性与可及性',
                'example': '(1)在社区卫生服务中心中(l),为提高2型糖尿病患者(c)的慢病管理质量(e)，由全科医生与护士团队(p)实施多学科管理模式(se)，以改善血糖控制率并加强并发症监测(i)',
                'elements': [
                    {'label': 'E', 'order': 1, 'description': 'Exceptation(improvement, innovation or information) - 期望'},
                    {'label': 'C', 'order': 2, 'description': 'Client group (recipients of service) - 服务对象'},
                    {'label': 'L', 'order': 3, 'description': 'Location (where service is housed) - 地点'},
                    {'label': 'I', 'order': 4, 'description': 'Impact (change in service and how measured) - 影响'},
                    {'label': 'P', 'order': 5, 'description': 'Professionals - 专业人员'},
                    {'label': 'Se', 'order': 6, 'description': 'Services - 服务'}
                ]
            }
        ]
        with self.connection.cursor() as cur:
            for framework in frameworks:
                cur.execute(insert_framework_query, (framework['label'], framework['supported_research_type'], framework['description']))
                result = cur.fetchone()
                if result:
                    framework_id = result[0]
                else:
                    cur.execute(search_framework_query, (framework['label'],))
                    framework_id = cur.fetchone()[0]

                for element in framework['elements']:
                    cur.execute(insert_framework_element_query, (framework_id, element['label'], element['order'], element['description']))  
        self.commit()

    def get_all_frameworks(self) -> List[Dict]:
        search_query = "SELECT * FROM framework ORDER BY id"
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query)
            result = cur.fetchall()
        self.commit()
        return [dict(row) for row in result] if result else []

    def get_framework_elements(self, framework_id: int) -> List[Dict]:
        search_query = """
                    SELECT * FROM framework_elements 
                    WHERE framework_id = %s
                    ORDER BY "order"
                    """
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (framework_id,))
            result = cur.fetchall()
        self.commit()
        return [dict(row) for row in result] if result else []

    def get_framework_element_by_label(self, framework_label: str, element_label: str) -> Optional[Dict]:
        search_query = """
                    SELECT fe.* 
                    FROM framework_elements fe
                    JOIN framework f ON fe.framework_id = f.id
                    WHERE f.label = %s AND fe.label = %s
                    """
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (framework_label, element_label))
            row = cur.fetchone()
        self.commit()
        return dict(row) if row else None

    def get_framework_by_label(self, label: str) -> Optional[Dict]:
        search_query = "SELECT * FROM framework WHERE label = %s"
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (label,))
            row = cur.fetchone()
        self.commit()
        return dict(row) if row else None
    
    def insert_question_source(self, source_label: str, research_type: str, grade: str, title: str, abstract: str, source_id: Optional[str] = None) -> Optional[int]:
        insert_query = """
                    INSERT INTO question_source 
                    (source_label, source_id, research_type, grade, title, abstract)
                    VALUES (%s, %s, %s, %s, %s, %s)
                    ON CONFLICT (source_label, title) DO NOTHING
                    RETURNING id
                    """
        with self.connection.cursor() as cur:
            cur.execute(insert_query, (source_label, source_id, research_type, grade, title, abstract))
            result = cur.fetchone()
        self.commit()
        return result[0] if result else None

    def get_question_source_by_id(self, question_source_id: int) -> Optional[Dict]:
        search_query = "SELECT * FROM question_source WHERE id = %s"
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (question_source_id,))
            row = cur.fetchone()
        self.commit()
        return dict(row) if row else None

    def get_question_source_by_title(self, source_label: str, title: str) -> Optional[Dict]:
        search_query = "SELECT * FROM question_source WHERE source_label = %s AND title = %s"
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (source_label, title))
            row = cur.fetchone()
        self.commit()
        return dict(row) if row else None
    
    def insert_question(self, question_source_id: int, framework_id: int, embedding: Optional[List[float]] = None) -> Optional[int]:
        insert_query = """
                INSERT INTO question 
                (question_source_id, framework_id, embedding)
                VALUES (%s, %s, null)
                ON CONFLICT (question_source_id, framework_id) DO NOTHING
                RETURNING id
                """
        with self.connection.cursor() as cur:
            cur.execute(insert_query, (question_source_id, framework_id))
            result = cur.fetchone()[0]
        if embedding is not None:
            update_query = """
                UPDATE question 
                SET embedding = %s
                WHERE id = %s
                """
            with self.connection.cursor() as cur:
                cur.execute(update_query, (embedding, result))
        self.commit()
        return result if result else None

    def get_question_by_id(self, question_source_id: int) -> Optional[List[Dict]]:
        search_query = "SELECT * FROM question WHERE question_source_id = %s"
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (question_source_id,))
            result = cur.fetchall()
        self.commit()
        return [dict(row) for row in result] if result else None

    def insert_question_element(self, question_id: int, framework_element_id: int, question_element_label: str) -> None:
        insert_query = """
                    INSERT INTO question_elements 
                    (question_id, framework_element_id, question_element_label)
                    VALUES (%s, %s, %s)
                    ON CONFLICT (question_id, framework_element_id) 
                    DO NOTHING
                    RETURNING id
                    """
        with self.connection.cursor() as cur:
            cur.execute(insert_query, (question_id, framework_element_id, question_element_label))
            result = cur.fetchone()
        self.commit()
        return result[0] if result else None

    def get_question_elements_by_question_id(self, question_id: int) -> Optional[List[Dict]]:
        search_query = """
                    SELECT 
                        qe.*,
                        fe.label as element_label,
                        fe.description as element_description
                    FROM question_elements qe
                    JOIN framework_elements fe ON qe.framework_element_id = fe.id
                    WHERE qe.question_id = %s
                    ORDER BY fe."order"
                    """
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (question_id,))
            result = cur.fetchall()
        self.commit()
        return [dict(row) for row in result] if result else None

    def get_statistics(self) -> Dict:
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            stats = {}
            
            cur.execute("SELECT COUNT(*) as count FROM framework")
            stats['framework_count'] = cur.fetchone()['count']
            
            cur.execute("SELECT COUNT(*) as count FROM question_source")
            stats['question_source_count'] = cur.fetchone()['count']
            
            cur.execute("SELECT COUNT(*) as count FROM question")
            stats['question_count'] = cur.fetchone()['count']
            
            cur.execute(
                """
                SELECT research_type, COUNT(*) as count 
                FROM question_source 
                GROUP BY research_type
                """
            )
            stats['by_research_type'] = {row['research_type']: row['count'] for row in cur.fetchall()}
            
            cur.execute(
                """
                SELECT f.label, COUNT(q.id) as count 
                FROM framework f
                LEFT JOIN question q ON f.id = q.framework_id
                GROUP BY f.label
                """
            )
            stats['by_framework'] = {row['label']: row['count'] for row in cur.fetchall()}
            
        self.commit()
        return stats
