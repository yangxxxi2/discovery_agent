from storage.db import DB
from dotenv import load_dotenv
from .knowledge_base_constructor import KnowledgeBaseConstructor

load_dotenv()

if __name__ == "__main__":
    
    db = DB()
    db.reset_db()

    knowledge_base = KnowledgeBaseConstructor()
    abstracts = knowledge_base.fetch_abstracts_from_pubmed(max_results=10)

    success_count = 0
    for item in abstracts:
        if(knowledge_base.persist_parse_abstract(item, db)):
            success_count += 1
    print(f"成功处理{success_count}篇摘要")

    db.close()
