from fastapi import FastAPI, HTTPException, Query
from typing import List, Optional
from .database import db
from .models import SampleInfo, ReadItem

app = FastAPI(
    title="Worf Genomics API",
    description="API for accessing WGS mapping data stored in MongoDB",
    version="1.0.0"
)

@app.get("/", tags=["General"])
def read_root():
    return {"status": "online", "system": "Worf Pipeline Server"}

@app.get("/samples", response_model=List[SampleInfo], tags=["Samples"])
def get_samples():
    """获取所有已入库的样品列表"""
    return db.get_samples()

@app.get("/query/{sample_id}", response_model=List[ReadItem], tags=["Query"])
def query_reads(
    sample_id: str,
    chrom: str,
    start: Optional[int] = Query(None, description="Start position (0-based)"),
    end: Optional[int] = Query(None, description="End position (exclusive)"),
    # 【修改点】将 le=5000 改为 le=100000 (十万) 或者直接去掉 le 限制
    limit: int = Query(100, ge=1, le=100000, description="Max records to return")
):
    """
    核心查询接口：根据染色体位置提取 Reads
    - 如果不填 start/end，则查询整条染色体（受 limit 限制）
    """
    results = db.query_reads(sample_id, chrom, start, end, limit)
    
    if not results:
        # 可以返回空列表，或者 404，这里选择空列表更友好
        return []
    
    return results