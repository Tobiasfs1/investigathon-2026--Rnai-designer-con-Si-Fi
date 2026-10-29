# Integración real con si-Fi y bases instaladas

Inspección de fuentes originales el 27/10/2026; inventario local por `rnai-designer databases --config examples/config.windows.json --inspect` y herramientas por `rnai-designer doctor`. Las bases y ejecutables originales se trataron solo en lectura. No se descargó, copió, reindexó ni modificó ningún índice.

## Qué se puede automatizar

**RNAi Designer: manual mediante la GUI.** `main.py:main` inicia QApplication. `MyMainWindow.start_sifi_pipeline` construye `SifiPipeline` y lee `run_pipeline`; `SifiPipeline.run_pipeline` abre la selección de main targets (`get_main_target` → `popup.ListSelection`) y luego abre `show_plot.DrawPlot`. No encontré CLI de si-Fi, biblioteca headless documentada ni forma estable de suprimir la GUI. Desde Python 3 sí se puede preparar el FASTA para abrir y ejecutar en RNAi Designer, con `rnai-designer prepare --candidate ... --out ...`; el alias de FASTA permite asociar su query ID al resultado. Luego `collect` recupera un único JSON temporal inequívocamente compatible o informa ambigüedad. No inventa métricas si si-Fi no se ejecutó.

**Búsquedas off-target: automatizadas con el Bowtie incluido en la instalación si-Fi.** El modo `Off-target Predictions` original es igualmente parte de `SifiPipeline` y por tanto tiene GUI. Para automatizar exactamente la búsqueda de secuencias, el módulo usa el mismo ejecutable Bowtie instalado y sus índices existentes, como proceso separado Python 3. La invocación probada es equivalente a `bowtie -a -v 0 -y <index> -f queries.fasta alignments.tsv`: todas las ventanas 21nt, cero mismatches, sin gaps y alineamientos de ambas hebras. No se importa Python 2 ni se reescribe un algoritmo de alineamiento. La automatización Bowtie **no ejecuta RNAi Designer ni RNAplfold**, y por eso sus registros no incluyen accesibilidad ni predicciones de eficacia.

El ensayo de integración ejecuta Bowtie real 0.11.3 contra índices pequeños creados bajo directorios temporales de prueba; prueba coincidencias exactas, hebra inversa, cero hits, reimportación y que los índices de prueba no cambian. No se ejecutó Bowtie contra las bases instaladas con candidatos biológicos: los únicos candidatos presentes son fixtures sintéticos. Para consultar bases reales, usar una secuencia e identificación aprobadas por el equipo.

## Instalación encontrada

La ruta se resuelve mediante `%LOCALAPPDATA%\siFi2015\Databases`; no se codifica el usuario. Todas las bases activas indicadas en examples/config.windows.json están en raíz con prefijos que coinciden exactamente con archivos observados:

| Organismo configurado | Prefijo | Formato | Estado |
|---|---|---|---|
| *Dalbulus maidis* | `Dalbulus_maidis` | Bowtie 1, seis `.ebwt` | Completo |
| *Homo sapiens* | `Homo_sapiens` | Bowtie 1, seis `.ebwt` | Completo |
| *Zea mays* | `zea_mays` | Bowtie 1, seis `.ebwt` | Completo |
| *Apis mellifera* | `Apis_mellifera` | Bowtie 1, seis `.ebwt` | Completo |

Los seis sufijos requeridos son `.1.ebwt`, `.2.ebwt`, `.3.ebwt`, `.4.ebwt`, `.rev.1.ebwt` y `.rev.2.ebwt`; tienen tamaño no cero. También hay una carpeta `Off_Targets/` con copias adicionales de índices de humano, maíz y abeja. Se enumera por separado; la configuración activa apunta a índices de raíz. No se eligió cuál copia adicional es preferible ni se sustituyeron las rutas configuradas.

`bowtie-inspect -n` tuvo código 0 para los cuatro índices raíz. Encontró 56.117 IDs en Dalbulus, 186.185 en Homo, 76.549 en Zea y 27.887 en Apis. IDs de muestra: Dalbulus `GITV01000001_1`; Homo, Zea y Apis devuelven IDs `NM_...` descritos como mRNA. Esto demuestra que los índices están legibles y da una pista de que varias referencias parecen transcritos; **no demuestra completitud taxonómica, ensamblado, accesiones de entrada ni fecha de la base**. No se encontró un manifiesto fiable de versión biológica. El código no declara que el contenido humano sea un genoma completo. Criterios negativos se limitan estrictamente a secuencias incluidas en el índice, no a todo el organismo.

`doctor` detectó:

| Ejecutable | Ruta resuelta desde LOCALAPPDATA | Versión |
|---|---|---|
| Bowtie | `siFi2015/Bowtie/bowtie.exe` | 0.11.3 (32-bit) |
| Bowtie inspector | `siFi2015/Bowtie/bowtie-inspect.exe` | ayuda disponible |
| RNAplfold | `siFi2015/RNAplfold/RNAplfold.exe` | 2.1.5 |

Estas versiones son antiguas, pero el Bowtie local ejecutó las pruebas contra índices sintéticos. `--config`, `--db-dir`, `SIFI_DB_DIR`, `SIFI_BOWTIE`, `SIFI_BOWTIE_INSPECT` y `SIFI_RNAPLFOLD` permiten portabilidad. La configuración no requiere modificar carpetas del usuario. El inventario calcula nombres, sufijos, tamaños, marcas de tiempo y completitud; el fingerprint en cada corrida calcula SHA-256 de cada archivo de índice para trazabilidad. No confundir fingerprint con identidad taxonómica.

## Criterio exacto de 21 nt

El comando `analyze` es explícito y consulta las tres especies; `--db-dir` permite indicar otro directorio. Se puede probar solo herramientas y bases con `doctor` y `databases --inspect`. La salida de análisis debe ir a un directorio nuevo, ajeno a los índices y al código original. Cada corrida conserva query, comando, timestamps, logs stdout/stderr, código/estado del proceso, conteo de ventanas, versión y hashes.

Una coincidencia se marca exacta solo si la fila original de Bowtie identifica la ventana de la query, la secuencia coincide con esa ventana (reverse complement para hebra `-`), el ID de referencia y offset son válidos, no hay descriptor de mismatch y no hay CIGAR que sugiera gaps. Los estados por organismo y global son:

- **Coincidencia exacta detectada:** existe al menos una fila verificable que contiene un 21-mer exacto.
- **Criterio computacional cumplido:** la ejecución completó, el archivo y la query verifican sus hashes, Bowtie buscó todos los 21-mers a cero mismatches, sin gaps, ambas hebras y todas las alineaciones, procesó todas las ventanas, y las filas importadas se interpretan correctamente sin coincidencias exactas.
- **Indeterminado:** los resultados no permiten demostrar ni coincidencia exacta ni búsqueda negativa completa, por ejemplo por parámetros ausentes, filas de mismatches no verificables, query ambigua o cobertura incompleta.
- **No evaluado:** la especie no tuvo análisis.
- **Error:** el proceso terminó mal, faltaron índices o el archivo quedó inválido.

Un alineamiento con mismatch se conserva y se distingue de un exacto; si la cobertura no es completa, los otros resultados permanecen indeterminados. El mensaje “criterio cumplido” se limita a los IDs de secuencia representados en esa base y los parámetros registrados. No significa inocuidad, no mide organismos/tejidos no representados y no excluye otras rutas de off-target. El uso real necesita que el equipo registre accesión y release de FASTA fuente de cada índice.

## RNAi Designer, JSON y errores

El modo original usa RNAplfold, strand selection, end stability, accesibilidad y reglas terminales documentadas en `SIFI_ANALYSIS.md`. Requiere tamaño de siRNA, base Bowtie, mismatches y configuración de filtros/umbrales. La selección manual de main targets es propia del modo RNAi design. En ambos modos la implementación original pasa por GUI y llama Bowtie/RNAplfold con procesos, pero abre selección/plots; no se verificó ejecución headless estable.

El JSON temporal es una lista generada por `SifiPipeline.data_to_json`; no guarda parámetros, organismo o versión de base. La copia inspeccionada filtra hits de hebra `-` y tiene un fallo al producir salida sin hits. Para datos originales del modo Off-target Predictions, los hits no se pueden reinterpretar como búsqueda negativa completa: el reporte los clasifica indeterminado por el formato/exportación. El adaptador de Bowtie conserva filas originales y es necesario para acreditar el criterio de ambas hebras. No se convierte un CSV resumen en evidencia posicional.

El modo RNAi sigue requiriendo: abrir el FASTA preparado, seleccionar explícitamente la base Dalbulus, seleccionar main targets, ejecutar si-Fi, completar `run.json` con parámetros y recoger el JSON correcto. `collect` no decide entre varios JSON compatibles; pide elegir uno explícitamente. Los off-target sí se lanzan directamente, con timeout por especie. Un error/timeout produce estado `error`, nunca cero hits ni criterio cumplido. Las métricas de eficiencia son las calculadas por si-Fi; el módulo no implementa optimización.
