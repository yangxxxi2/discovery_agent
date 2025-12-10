from .constants import TASK_TO_MODEL_MAP
from sentence_transformers import SentenceTransformer
import re
from .llm import llm
import json
from storage.db import DB
from .prompt import DETERMINE_TASK_TYPE_PROMPT, EXTRACT_QUESTION_MODEL_PROMPT, EXTRACT_RESEARCH_TYPE_PROMPT
from .utils import format_components_info, format_example_json
import numpy as np

class QuestionModelExtractor:

    def determine_task_type(self, abstract: dict) -> str:
        try:
            title = abstract.get("title", "")
            abstract_text = abstract.get("abstract", "")
            formatted_prompt = DETERMINE_TASK_TYPE_PROMPT.format(title=title, abstract=abstract_text)
            
            result = llm.invoke(formatted_prompt)
            result_text = result.content if hasattr(result, 'content') else str(result)
            valid_types = list(TASK_TO_MODEL_MAP.keys())
            if result_text in valid_types:
                return result_text
            for valid_type in valid_types:
                if valid_type.lower() in result_text.lower() or result_text.lower() in valid_type.lower():
                    return valid_type
            return None
        
        except Exception as e:
            print(f"判断任务类型时出错: {e}")
            return None
        

    def extract_question_model(self, abstract: dict, task_type: str, db: DB) -> dict:
        model_type = TASK_TO_MODEL_MAP.get(task_type, None)
        if not model_type:
            return {}
        try:
            model_id = db.get_framework_by_label(model_type).get("id", None)
            model = db.get_framework_by_id(model_id)
            if not model_id:
                print(f"未找到模型 {model_type} ")
                return {}
            components_list = db.get_framework_elements(model_id)
            components_info = format_components_info(components_list)
            example_format = format_example_json(components_list)

            title = abstract.get("title", "")
            abstract = abstract.get("abstract", "")
            formatted_prompt = EXTRACT_QUESTION_MODEL_PROMPT.format(
                model_type=model_type, 
                title=title, 
                abstract=abstract, 
                components_info=components_info, 
                example_format=example_format, 
                model_description=model.get("description", ""), 
                model_examples=model.get("example", ""))
            
            result = llm.invoke(formatted_prompt)
            result = result.content if hasattr(result, 'content') else json.loads(str(result))
            return {"model_type": model_type, **result}
        
        except json.JSONDecodeError as e:
            print(f"解析LLM返回的JSON失败: {e}")
            return None
        except Exception as e:
            print(f"抽取问题模型时出错: {e}")
            return {}


    def extract_research_type(self, abstract: dict) -> str:
        try:
            title = abstract.get("title", "")
            abstract_text = abstract.get("abstract", "")
            formatted_prompt = EXTRACT_RESEARCH_TYPE_PROMPT.format(title=title, abstract=abstract_text)

            result = llm.invoke(formatted_prompt)
            result_text = result.content if hasattr(result, 'content') else str(result)
            result_text = result_text.strip()
            
            valid_research_types = [
                'Randomized Controlled Trial',
                'Non Randomized Controlled Trial',
                'Single Arm Trial',
                'Three Arm Trial',
                'Crossover Trial',
                'Factorial Trial',
                'Cohort Study',
                'Case-Control Study',
                'Cross-Sectional Study',
                'Before-and-After Study',
                'Ecological Study',
                'Focus Group Discussion',
                'In-depth Interview'
            ]
            if result_text in valid_research_types:
                return result_text
            for valid_type in valid_research_types:
                if valid_type.lower() in result_text.lower() or result_text.lower() in valid_type.lower():
                    return valid_type
                return ""
        except Exception as e:
            print(f"提取研究设计类型时出错: {e}")
            return ""

    def embed_question_model(self, question_model: dict) -> np.ndarray:
        try:
            texts = []
            for key, value in question_model.items():
                if value != "Not specified" and value.strip():
                    texts.append(f"{key}: {value}")
            combined_text = re.sub(r"[\n\r\t]", "", "\n".join(texts)).strip()
            full_text = f"Model: {question_model.model_type}; {combined_text}"
            embedding_model = SentenceTransformer('BAAI/bge-m3', device='cpu')
            vector = embedding_model.encode([full_text], convert_to_tensor=False, normalize_embeddings=True)
            return vector
        except Exception as e:
            print(f"生成embedding时出错: {e}")
            return np.array([])
