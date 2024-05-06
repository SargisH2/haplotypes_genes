import uvicorn
from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from gene_haplotype_functions import get_genes, get_hla

app = FastAPI()

@app.get("/")
async def root():
    return {"message": 'Send post query: {"text": "<<text>>"}'}

@app.post('/')
async def get_haplotypes(request: Request):
    data = await request.json()
    text = data.get('text')
    if not isinstance(text, str):
        return JSONResponse(status_code=400, content={"error": "Input should be in text format"})
    
    return JSONResponse(
        status_code=200,
        content = {
            "hla": get_hla(text), 
            "genes": get_genes(text)
        }
    )

if __name__ == "__main__":
    uvicorn.run("web_service:app", host="0.0.0.0", port=8070)
