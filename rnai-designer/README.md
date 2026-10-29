# rnai-designer

Capa Python 3 de automatización e integración con si-Fi para investigación de RNAi en *Dalbulus maidis*. Versión declarada: **0.3.0**. Conserva las métricas originales y añade ventanas, prioridad heurística, selección condicionada al filtro off-target, batch y reportes trazables.

**La prioridad no es eficacia experimental ni garantía de seguridad biológica.**

## Guía del proyecto

- [README principal](../README.md): presentación, pipeline, carpetas y limitaciones.
- [Índice técnico](../docs/README.md).
- [Instalación y versiones verificadas](../docs/instalacion-y-reproducibilidad.md).
- [Arquitectura](../docs/arquitectura.md).
- [Métricas, puntuación y off-target](../docs/criterios-cientificos.md).
- [Batch, ranking y reportes](../docs/batch-y-reportes.md).
- [Rendimiento y crecimiento de archivos](../docs/rendimiento-y-limitaciones.md).
- [Créditos y licencias](../docs/creditos-y-licencias.md).

## Inicio rápido en PowerShell

Desde esta carpeta, con Python 3.10 o posterior:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m rnai_designer doctor
```

El workflow requiere además el entorno Conda Python 2.7 de si-Fi, Bowtie clásico, RNAplfold y los índices. La instalación pip no prepara esos recursos. No instalar ni ejecutar diagnósticos sobre un entorno que tenga un batch activo.

## Ejecución

Con el Python 3 configurado, o sustituyendo `python` por `.\.venv\Scripts\python.exe`:

```powershell
python -m rnai_designer MiGen "SECUENCIA_IUPAC_COMPLETA"
python -m rnai_designer batch genes.fasta
python -m rnai_designer batch --resume output/batch_1
```

Reemplazar el marcador por nucleótidos reales. El comando individual crea FASTA y carpeta corta correlativa sin sobrescribirlos. El batch conserva entradas duplicadas por posición, continúa ante fallos individuales y genera `resultados_genes.csv`, `ranking_genes.csv` y `resumen_batch.json`.

En batch, cada gen tiene un límite de **360 segundos (6 minutos)**. Usar `--gene-timeout 600` para cambiarlo o `--gene-timeout 0` para desactivarlo. Los vencidos quedan como `incomplete_timeout`, sin puntaje en el ranking, y pueden reintentarse con `batch --resume`. La cancelación incluye los procesos descendientes del trabajador. El comando individual `workflow` conserva sus mecanismos existentes.

Forma avanzada, con archivo existente y salida nueva:

```powershell
python -m rnai_designer workflow --gene-id MiGen --sequence-file inputs/MiGen.fasta --out output/MiGen_revision --run-sifi-headless
```

El resumen es compacto por defecto; `workflow --verbose-report` habilita el detalle. El puntaje original pondera 75 % eficiencia si-Fi y 25 % GC. La diversidad conserva su cálculo, exportación y papel en los desempates, con 0 % de peso en el puntaje. Las ventanas se prueban por puntaje hasta acreditar el filtro en las tres especies. Si ninguna cumple, no hay ganadora ni secuencia publicada. Los puntajes originales no se penalizan por off-target; su fórmula y la penalización por longitud se documentan por separado.

## Documentos especializados previos

Se conservan [integración](docs/SIFI_INTEGRATION.md), [formatos](docs/FORMATS.md), [auditoría original](docs/SIFI_ANALYSIS.md), [automatización](docs/SIFI_AUTOMATION.md), [portabilidad](docs/portabilidad.md), [ventanas](docs/WINDOW_EXPLORATION.md), [selección](docs/AUTOMATIC_WINDOW_SELECTION.md), [optimización](docs/PERFORMANCE_OPTIMIZATION.md) y [validación](docs/VALIDATION.md).

Pueden describir etapas anteriores del desarrollo. Para el comportamiento actual, utilizar la guía técnica del proyecto y el código correspondiente al commit analizado. Los datos de `examples/` son ficticios y no constituyen evidencia biológica.
