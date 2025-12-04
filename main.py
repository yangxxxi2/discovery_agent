from storage.db import DB
from dotenv import load_dotenv

load_dotenv()

if __name__ == "__main__":
    
    db = DB()
    db.reset_db()
    
    stats = db.get_statistics()
    print("\n数据库统计信息:")
    for key, value in stats.items():
        print(f"{key}: {value}")