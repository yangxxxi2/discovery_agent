from typing import List, Dict
from .question_model_extractor import QuestionModelExtractor
from Bio import Entrez
from .constants import KEYWORD_CLINICAL_RESEARCH, TASK_TO_MODEL_MAP
from xml.etree import ElementTree as ET
import re
import time
from storage.db import DB

class KnowledgeBaseConstructor:

    def __init__(self):
        self.extractor = QuestionModelExtractor()

    def fetch_abstracts_from_pubmed(self, max_results: int) -> List[Dict]:
        try:
            abstracts = []

            handle = Entrez.esearch(db="pubmed", term=KEYWORD_CLINICAL_RESEARCH, retmax=max_results)
            record = Entrez.read(handle)
            pmids = record["IdList"]
            time.sleep(0.5)

            if not pmids:
                print(f"未找到关键词 '{KEYWORD_CLINICAL_RESEARCH}' 相关的文献")
                return []

            for pmid in pmids:
                time.sleep(0.5)
                try:
                    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="xml")
                    xml_data = handle.read()
                    root = ET.fromstring(xml_data)

                    title = root.find(".//ArticleTitle").text
                    if not title:
                        continue

                    for abstract in root.findall(".//Abstract"):
                        for abstract_text_elem in abstract.findall("AbstractText"):
                            text += abstract_text_elem.text
                    clean_text = re.sub(r"[\n\r\t]", "", text).strip()
                    if not clean_text or len(clean_text) < 100:
                        continue

                    abstracts.append({"pmid": pmid, "abstract": clean_text, "title": title})
            
                except Exception as e:
                    print(f"解析文献PMID {pmid} 的摘要和标题时出错: {e}")
                    continue
        
            print(f"成功获取 {len(abstracts)} 篇文献摘要")
            return abstracts
            
        except Exception as e:
            print(f"从pubmed获取摘要时发生错误: {e}")
            return []

    def _parse_abstract(self, abstract: Dict, db: DB) -> Dict:
        title = abstract.get("title", "")
        abstract_text = abstract.get("abstract", "")
        text = "title: " + title + "\nabstract: " + abstract_text
        task_type = self.extractor.determine_task_type(text)
        question_model = self.extractor.extract_question_model(text, task_type, db) if task_type != None else {}
        research_type = self.extractor.extract_research_type(text) if task_type != None else ""
        embedded_model = self.extractor.embed_question_model(question_model) if question_model != {} else []
        parsed_abstract = {
                **question_model,
                "task_type": task_type,
                "research_type": research_type,
                "embedding": embedded_model,
            }
        return parsed_abstract

    def persist_parse_abstract(self, item: dict, db: DB) -> bool:
        parsed_abstract = self._parse_abstract(item)
        if len(parsed_abstract["embedding"]) == 0 | parsed_abstract["question_model"] == {} | parsed_abstract["research_type"] == "":
            print("解析失败，跳过存储")
            return False
        try:
            question_source_id = db.insert_question_source(
                source_label="PubMed",
                source_id=item["pmid"],
                research_type=parsed_abstract["research_type"],
                task_type=parsed_abstract["task_type"],
                grade= "",
                ocebm_2011_level= "",
                title=item["title"],
                abstract=item["abstract"]
            )
            if question_source_id == None:
                return False
            model_type = TASK_TO_MODEL_MAP.get(parsed_abstract["task_type"], None)
            framework_id = db.get_framework_by_label(model_type).get("id", None)
            if framework_id == None:
                return False
            model_id = db.insert_question(
                question_source_id=question_source_id,
                framework_id=framework_id,
                embedding=parsed_abstract["embedding"],
            )
            if model_id == None:
                return False
            for key, value in parsed_abstract.items():
                if key not in ["task_type", "research_type", "embedding"]:
                    element_id = db.get_framework_element_by_label(framework_id, key).get("id", None)
                    if element_id:
                        db.insert_question_component(
                            question_id=model_id,
                            framework_element_id=element_id,
                            question_element_label=value
                        )
                        db.commit()
                        return True
                    else:
                        return False
                else:
                    return False
        except Exception:
            print(f"持久化模型到数据库中出错")
            return False
