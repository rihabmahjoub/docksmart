"""DockSmart — FastAPI application entrypoint."""
from __future__ import annotations

import logging

from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from app.config import settings
from app.database import init_db
from app.routers import jobs, structures

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

app = FastAPI(
    title=settings.APP_NAME,
    version=settings.APP_VERSION,
    description=(
        "DockSmart: an automated web tool for ligand-protein interaction studies. "
        "Automates receptor/ligand preparation, binding-site detection, AutoDock "
        "Vina docking, and result interpretation."
    ),
)

app.mount("/static", StaticFiles(directory="app/static"), name="static")
templates = Jinja2Templates(directory="app/templates")

app.include_router(jobs.router)
app.include_router(structures.router)


@app.on_event("startup")
def on_startup() -> None:
    init_db()


@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    return templates.TemplateResponse(request, "index.html")


@app.get("/jobs/{job_id}", response_class=HTMLResponse)
def job_page(request: Request, job_id: str):
    return templates.TemplateResponse(request, "job_status.html", {"job_id": job_id})


@app.get("/health")
def health():
    return {"status": "ok", "app": settings.APP_NAME, "version": settings.APP_VERSION}
