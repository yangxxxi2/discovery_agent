import os
import psycopg2
from psycopg2.extras import RealDictCursor
import logging
from typing import List, Dict, Optional

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

    def create_tables(self, enable_vector: bool = True):
        """创建数据库表
        
        Args:
            enable_vector: 是否启用pgvector扩展。如果为False，embedding列将使用TEXT类型
        """
        # 尝试创建vector扩展，如果失败则使用TEXT类型
        vector_type = "VECTOR(1536)" if enable_vector else "TEXT"
        
        if enable_vector:
            try:
                with self.connection.cursor() as cur:
                    cur.execute("CREATE EXTENSION IF NOT EXISTS vector;")
                self.connection.commit()
                logger.info("pgvector扩展已启用")
            except Exception as e:
                logger.warning(f"pgvector扩展启用失败，将使用TEXT类型: {e}")
                vector_type = "TEXT"
                enable_vector = False
                self.connection.rollback()
        
        create_query = f"""

        CREATE TABLE IF NOT EXISTS framework (
            id                      SERIAL PRIMARY KEY,
            label                   VARCHAR(50) UNIQUE NOT NULL,
            supported_research_type VARCHAR(100) NOT NULL,
            description             TEXT,
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
            grade                   VARCHAR(50) NOT NULL,
            title                   TEXT NOT NULL,
            abstract                TEXT NOT NULL,
            CONSTRAINT unique_source UNIQUE(source_label, title)
        );

        CREATE TABLE IF NOT EXISTS question (
            id                      SERIAL PRIMARY KEY,
            question_source_id      INT NOT NULL REFERENCES question_source(id) ON DELETE CASCADE,
            framework_id            INT NOT NULL REFERENCES framework(id),
            embedding               {vector_type} NOT NULL,
            metadata                JSONB,
            CONSTRAINT unique_question UNIQUE(question_source_id, framework_id)
        );

        CREATE TABLE IF NOT EXISTS question_elements (
            id                      SERIAL PRIMARY KEY,
            question_id             INT NOT NULL REFERENCES question(id) ON DELETE CASCADE,
            framework_element_id    INT NOT NULL REFERENCES framework_elements(id),
            question_element_label  TEXT NOT NULL,
            CONSTRAINT unique_question_element UNIQUE(question_id, framework_element_id)
        );

        """
        
        with self.connection.cursor() as cur:
            cur.execute(create_query)
        
        # 只有在启用vector时才创建向量索引
        if enable_vector:
            try:
                with self.connection.cursor() as cur:
                    cur.execute("""
                        CREATE INDEX IF NOT EXISTS idx_question_embedding_vector 
                        ON question USING ivfflat (embedding vector_cosine_ops) 
                        WITH (lists = 100);
                    """)
                logger.info("向量索引创建成功")
            except Exception as e:
                logger.warning(f"向量索引创建失败: {e}")
                self.connection.rollback()
        
        self.connection.commit()
        logger.info("数据表创建成功")

    def drop_tables(self):
        """删除所有表"""
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
        self.connection.commit()
        logger.info("所有数据表已删除")

    def initialize_default_frameworks(self) -> None:
        insert_framework_query = """
                        INSERT INTO framework (label, supported_research_type, description)
                        VALUES (%s, %s, %s)
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
                'label': 'PICO',
                'supported_research_type': 'Randomized Controlled Trials',
                'description': 'Patient/Population, Intervention, Comparison, Outcome',
                'elements': [
                    {'label': 'P', 'order': 1, 'description': 'Population - 研究人群'},
                    {'label': 'I', 'order': 2, 'description': 'Intervention - 干预措施'},
                    {'label': 'C', 'order': 3, 'description': 'Comparison - 对照组'},
                    {'label': 'O', 'order': 4, 'description': 'Outcome - 结局指标'}
                ]
            },
            {
                'label': 'SPIDER',
                'supported_research_type': 'Qualitative Research',
                'description': 'Sample, Phenomenon of Interest, Design, Evaluation, Research type',
                'elements': [
                    {'label': 'S', 'order': 1, 'description': 'Sample - 样本/研究对象'},
                    {'label': 'PI', 'order': 2, 'description': 'Phenomenon of Interest - 关注现象'},
                    {'label': 'D', 'order': 3, 'description': 'Design - 研究设计'},
                    {'label': 'E', 'order': 4, 'description': 'Evaluation - 评估/结果'},
                    {'label': 'R', 'order': 5, 'description': 'Research type - 研究类型'}
                ]
            },
            # {
            #     'label': 'PICOC',
            #     'supported_research_type': 'Software Engineering',
            #     'description': 'Population, Intervention, Comparison, Outcome, Context',
            #     'elements': [
            #         {'label': 'P', 'order': 1, 'description': 'Population - 研究对象'},
            #         {'label': 'I', 'order': 2, 'description': 'Intervention - 干预措施'},
            #         {'label': 'C', 'order': 3, 'description': 'Comparison - 对照'},
            #         {'label': 'O', 'order': 4, 'description': 'Outcome - 结果'},
            #         {'label': 'C', 'order': 5, 'description': 'Context - 上下文'}
            #     ]
            # },
            # {
            #     'label': 'SPICE',
            #     'supported_research_type': 'Mixed Methods Research',
            #     'description': 'Setting, Perspective, Intervention, Comparison, Evaluation',
            #     'elements': [
            #         {'label': 'S', 'order': 1, 'description': 'Setting - 研究背景'},
            #         {'label': 'P', 'order': 2, 'description': 'Perspective - 研究视角'},
            #         {'label': 'I', 'order': 3, 'description': 'Intervention - 干预'},
            #         {'label': 'C', 'order': 4, 'description': 'Comparison - 比较'},
            #         {'label': 'E', 'order': 5, 'description': 'Evaluation - 评估'}
            #     ]
            # }
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
        self.connection.commit()
        logger.info("默认框架初始化成功")

    def get_all_frameworks(self) -> List[Dict]:
        """获取所有框架"""
        search_query = "SELECT * FROM framework ORDER BY id"
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query)
            result = cur.fetchall()
        return [dict(row) for row in result] if result else []

    def get_framework_elements(self, framework_id: int) -> List[Dict]:
        """ 获取框架元素"""
        search_query = """
                    SELECT * FROM framework_elements 
                    WHERE framework_id = %s
                    ORDER BY "order"
                    """
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (framework_id,))
            result = cur.fetchall()
        return [dict(row) for row in result] if result else []

    def get_framework_element_by_label(self, framework_label: str, element_label: str) -> Optional[Dict]:
        """根据标签获取框架元素"""
        search_query = """
                    SELECT fe.* 
                    FROM framework_elements fe
                    JOIN framework f ON fe.framework_id = f.id
                    WHERE f.label = %s AND fe.label = %s
                    """
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (framework_label, element_label))
            row = cur.fetchone()
        return dict(row) if row else None

    def get_framework_by_label(self, label: str) -> Optional[Dict]:
        """根据标签获取框架"""
        search_query = "SELECT * FROM framework WHERE label = %s"
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (label,))
            row = cur.fetchone()
        return dict(row) if row else None
    
    def insert_question_source(self, source_label: str, research_type: str, grade: str, title: str, abstract: str, source_id: Optional[str] = None) -> Optional[int]:
        """插入问题来源"""
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
        self.connection.commit()
        return result[0] if result else None

    def get_question_source_by_id(self, question_source_id: int) -> Optional[Dict]:
        """根据ID获取问题来源"""
        search_query = "SELECT * FROM question_source WHERE id = %s"
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (question_source_id,))
            row = cur.fetchone()
        return dict(row) if row else None

    def get_question_source_by_title(self, source_label: str, title: str) -> Optional[Dict]:
        """根据标题获取问题来源"""
        search_query = "SELECT * FROM question_source WHERE source_label = %s AND title = %s"
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (source_label, title))
            row = cur.fetchone()
        return dict(row) if row else None
    
    def insert_question(self, question_source_id: int, framework_id: int, embedding: Optional[List[float]] = None) -> Optional[int]:
        """插入问题"""
        insert_query = """
                INSERT INTO question 
                (question_source_id, framework_id, embedding)
                VALUES (%s, %s, %s)
                ON CONFLICT (question_source_id, framework_id) DO NOTHING
                RETURNING id
                """
        # 如果embedding为list，转换为JSON字符串存储
        embedding_str = str(embedding) if embedding else None
        
        with self.connection.cursor() as cur:
            cur.execute(insert_query, (question_source_id, framework_id, embedding_str))
            result = cur.fetchone()
        self.connection.commit()
        return result[0] if result else None

    def get_question_by_id(self, question_source_id: int) -> Optional[List[Dict]]:
        """根据ID获取问题"""
        search_query = "SELECT * FROM question WHERE question_source_id = %s"
        with self.connection.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(search_query, (question_source_id,))
            result = cur.fetchall()
        return [dict(row) for row in result] if result else None

    def insert_question_element(self, question_id: int, framework_element_id: int, question_element_label: str) -> None:
        """插入问题元素"""
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
        self.connection.commit()
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
            
        return stats
