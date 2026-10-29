# Contratos de entrada y salida

UTF-8, JSON estricto sin NaN ni Infinity. Los ejemplos completos están en `examples/`. El candidato usa esquema 1.0 y el reporte integrado esquema 1.1. Los campos adicionales se conservan en candidato, registros y manifiestos.

## Candidato

Objeto con `schema_version: "1.0"`, `gene_id` (Dalbulus), `reference.sequence` (IUPAC), `metadata` (objeto) y `evidence` (lista). `tribolium_homolog_id`, `reference.source`, `reference.version` y `priority_score` pueden ser null/desconocidos al iniciar por secuencia; si están disponibles, se conservan. La prioridad se conserva sin recalcularla. Null significa desconocida, nunca cero implícito.

La referencia debe ser exactamente la query de la ejecución. Para un fragmento se debe suministrar su secuencia exacta y guardar el mapeo al gen en metadata. No se convierten coordenadas de fragmentos al gen completo. Se preserva la secuencia original; solo la comparación para advertencias ignora mayúsculas/minúsculas.

## Resultado si-Fi y manifiesto primario

Lista de registros verificada en SIFI_ANALYSIS.md. Se requieren `sirna_name`, `sirna_sequence`, `sirna_position` y coincidencia de `query_name` con gene_id o con query_id del manifiesto. Campos opcionales ausentes/null generan advertencias; tipos incompatibles, secuencias inválidas y coordenadas imposibles generan error. No se convierten strings booleanos ni el descriptor mismatches. Secuencias discordantes con la referencia se conservan con advertencia; deben revisarse antes de interpretar el gráfico. Null no se convierte a False.

`--run` acepta un objeto opcional con `query_id`, `tool`, `tool_version`, `source`, `database`, `database_version`, `parameters` y metadatos adicionales. Los campos de procedencia desconocidos pueden ser null u omitirse en el manifiesto primario y generan advertencias. Registrar tamaño, mismatches permitidos, filtros, umbrales, ventana, parámetros RNAplfold, fecha y accesión/versiones de referencia. No se atribuyen defaults a datos importados. Los hashes SHA-256 identifican contenido, no validan autoría ni calidad científica.

## Evaluaciones externas

`--offtargets` admite `schema_version: "1.0"` o `"1.1"` y `evaluations: [...]`. Se admite una evaluación por organismo con nombre exacto `Homo sapiens`, `Zea mays` o `Apis mellifera`. Cada entrada requiere `organism`, `query_id`, `tool`, `tool_version`, `database`, `database_version`, `source` (textos), `parameters` (objeto) y `results_file` (ruta relativa al manifiesto). El manifiesto 1.1 puede declarar `results_format: "bowtie_legacy_tab"`, ejecución, cobertura, versión y fingerprint de índices. La CLI `analyze` produce ese formato con Bowtie externo.

El formato 1.0 predeterminado es la lista JSON original de si-Fi. El formato `bowtie_legacy_tab` se importa con el adaptador de Bowtie, preserva las filas/tablas originales e interpreta coincidencias exactas por secuencia, offset y hebra. El CSV resumen de si-Fi y otros alineadores requieren adaptadores específicos: renombrar un archivo no convierte su formato.

Organismo omitido: `not_evaluated`. Archivo ausente: error. Archivo vacío: importado con advertencia salvo que un manifiesto completo acredite búsqueda exhaustiva. `imported_unvalidated` indica éxito técnico de importación, nunca validación de taxonomía de la base. Una evaluación negativa para 21nt exactos solo recibe `criterion_met` con formato Bowtie original, retorno 0, hashes concordantes de query/resultados, secuencia A/C/G/T, mismatches cero, todas las ventanas procesadas, búsqueda exhaustiva de alineamientos en ambas hebras y filas interpretables de 21nt. Mismatches se cuentan y registran aparte; filas dudosas dan `indeterminate`. Una coincidencia exacta verificada puede reportarse aun si la corrida quedó incompleta. Estados: `criterion_met`, `exact_match_detected`, `indeterminate`, `not_evaluated`, `error`. `safety_conclusion` siempre es `undetermined`.

## Salidas

- `report.json`: candidato y puntaje originales, registros, manifiestos, hashes, resumen descriptivo, tres organismos y advertencias.
- `records.csv`: una fila por hit con gene_id, prioridad, contexto y organismo; columnas métricas para análisis directo y `record_json` para preservar campos originales y extensiones. También conserva procedencia, parámetros y advertencias como celdas JSON.
- `evaluations.csv`: estado de ejecución y criterio de cada especie, con evaluación completa como celda JSON, incluso si no hubo filas.
- `offtarget_matches.csv`: ID de referencia, posición query y referencia, hebra, longitud, mismatches y coincidencia exacta verificada para cada alineamiento.
- `offtargets.svg` y `offtargets.png`: coincidencias importadas por posición y especie, con sus estados explícitos.
- La subcomanda `windows` crea en un directorio nuevo `window_analysis.json`, `windows.csv`, `window_offtargets.csv`, `windows.svg` y `windows.png`; `selected_window.json` se incluye si se especificó una ventana. Definición y campos: [WINDOW_EXPLORATION.md](WINDOW_EXPLORATION.md).
- `positions.csv`: puntos descriptivos deduplicados por siRNA e intervalo, con datos faltantes explícitos.
- `efficiency.svg` y `efficiency.png`: posición inicial 1-based contra accessibility_value y un segundo panel con is_efficient. No interpolan. En datos simulados se incluye aviso visible.

Los hits repetidos no inflan puntos. Métricas discordantes para una misma ventana producen advertencias y valores desconocidos en el resumen; el original permanece intacto. No se calcula un nuevo ranking biológico. La versión del contrato cambiará si cambia el significado o la obligatoriedad de campos; pueden agregarse metadatos preservando compatibilidad.
