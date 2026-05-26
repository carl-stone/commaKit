# Docker Container Package Requests

commaBot logs runtime package installs here. Carl reviews and either bakes them into the next image rebuild or explains why not.

## R packages

(none yet)

## System packages (apt)

- `pandoc` — needed for devtools::build_vignettes(). Currently missing from the container, which blocks vignette builds.
- `qpdf` — needed by `R CMD check` for PDF size-reduction checks when vignettes are built/compacted; without it local full checks report `qpdf is needed for checks on size reduction of PDFs`.

## Python packages

(none yet)

## Node.js packages

(none yet)

---

## How this works

1. commaBot installs a package at runtime (ephemeral — lost on container recreation)
2. commaBot adds it to the appropriate section above
3. Carl reviews, then either:
   - Adds it to the Dockerfile / renv.lock and rebuilds, then removes the entry
   - Tells commaBot why it won't be added and removes the entry
4. Next `docker compose up -d --build` picks up the change

## Rebuild reminder

```bash
# Sync renv files if R packages changed
cp /home/carls/comma/renv.lock /home/carls/commabot-infrastructure/comma/ && \
cp /home/carls/comma/renv/settings.json /home/carls/commabot-infrastructure/comma/renv/ && \
cp /home/carls/comma/renv/activate.R /home/carls/commabot-infrastructure/comma/renv/

# Rebuild
cd /home/carls/commabot-infrastructure && docker compose up -d --build
```
