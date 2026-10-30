# Arquitectura e integración con si-Fi

[Índice](README.md) · [Criterios científicos](criterios-cientificos.md) · [Rendimiento](rendimiento-y-limitaciones.md)

## Dos entornos separados

`rnai-designer` coordina el análisis en Python 3. El motor si-Fi se ejecuta en un subproceso Python 2.7 mediante `conda run -p`. No se importa Python 2 dentro del intérprete Python 3. El modo headless evita abrir la GUI, aunque las importaciones originales todavía necesitan componentes gráficos legados.

El motor está en [siFi21--1.2.3](../siFi21--1.2.3/README.md) y el adaptador en [legacy_bridge/run_sifi_headless.py](../rnai-designer/rnai_designer/legacy_bridge/run_sifi_headless.py). El adaptador llama métodos originales; las optimizaciones de integración no equivalen a una nueva validación biológica.

## Módulos principales

Todos los módulos siguientes pertenecen a `rnai-designer/rnai_designer/`:

| Módulo | Responsabilidad |
|---|---|
| `__main__.py` | CLI, coordinación por gen y main target |
| `simple_run.py`, `batch.py` | Nombres cortos, entradas múltiples y reanudación |
| `data_models.py`, `sifi_adapter.py` | Contratos, normalización, validación y procedencia |
| `config.py`, `legacy_runtime.py` | Herramientas, índices, entorno Conda y diagnóstico |
| `handoff.py`, `legacy_bridge/` | FASTA de intercambio y motor original |
| `efficiency_analysis.py`, `window_explorer.py` | Resúmenes y máscara verde ordenada |
| `automatic_window_selection.py`, `window_index.py`, `compact_windows.py` | Ventanas y agregación compacta de métricas |
| `priority_scoring.py` | Penalización por longitud y desempates |
| `automation.py`, `bowtie_adapter.py`, `offtarget_assessment.py` | Ejecución off-target, parsing y evaluación de evidencia |
| `safe_window_selection.py` | Recorrido del ranking hasta acreditar el filtro |
| `pipeline_workflow.py` | Reporte integrado y selección final |
| `reporting.py`, `automatic_window_reporting.py`, `console_output.py` | JSON, CSV, figuras y consola |
| `legacy_session.py`, `reuse.py`, `cache.py`, `atomic_io.py` | Cachés de distinto alcance y publicación atómica |

## Flujo por gen

1. Validar y normalizar la secuencia; conservar el nombre original y preparar el handoff con alias interno estable.
2. Detectar main targets mediante Bowtie contra la base de *D. maidis*.
3. Compartir los alineamientos, RNAplfold y registros base invariantes dentro de la ejecución del gen.
4. Para cada GITV, generar sus etiquetas `is_off_target`, JSON de si-Fi y máscara original; conservar el orden de filas.
5. Elaborar el reporte integrado y la evaluación off-target del gen completo.
6. Enumerar ventanas y calcular las métricas de [criterios científicos](criterios-cientificos.md).
7. Evaluar candidatas por prioridad hasta encontrar una que cumpla las tres especies o agotar las elegibles.
8. Exportar y consolidar. En batch se actualiza además el estado y el ranking general.

No se descartan GITV por pocos hits. Una referencia sin zonas verdes conserva sus reportes y no construye índices de ventanas ni evalúa una candidata inexistente.

## Bowtie y RNAplfold

**Bowtie** alinea secuencias cortas: encuentra referencias en la base primaria y coincidencias contra las bases off-target. No calcula accesibilidad ni eficacia. Se utiliza la familia clásica de Bowtie, no Bowtie 2.

**RNAplfold**, de ViennaRNA, aporta información de accesibilidad estructural utilizada por los criterios de eficiencia originales. No demuestra accesibilidad in vivo. Para la misma secuencia y parámetros, sus resultados se comparten entre GITV.

## Trazabilidad y formatos

Se conservan JSON originales, máscaras, parámetros, hashes, manifiestos y evidencias. El reporte integrado declara esquema JSON `1.1`; resúmenes utilizan `1.0` y tablas actuales CSV `2.0`. Consulte cada campo de esquema: no todas las salidas tienen la misma versión.

Las escrituras atómicas y el bloqueo del batch reducen inconsistencias, pero no hacen transaccional al conjunto entero de archivos. Una carpeta con archivos puede seguir siendo un intento parcial.
