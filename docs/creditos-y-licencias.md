# Créditos, componentes de terceros y licencias

[Índice](README.md) · [Presentación](../README.md) · [Herramientas](instalacion-y-reproducibilidad.md)

## si-Fi original

El [README original](../siFi21--1.2.3/README.md) identifica **siFi: Software for long double-stranded RNAi-target design and off-target prediction**, versión **v1.2.3-0008**, y acredita a:

- Stefanie Lück.
- Tino Kreszies.
- Marc Strickert.
- Patrick Schweizer.
- Dimitar Douchkov.

Ese mismo archivo declara **Attribution-NonCommercial-ShareAlike 2.0 Generic (CC BY-NC-SA 2.0)**. Esta documentación registra la declaración existente; no la sustituye por una licencia nueva ni determina por sí sola qué archivos derivados o recursos de terceros quedan cubiertos.

No se encontró un archivo de licencia independiente que ampliara esa declaración en la copia inspeccionada. La disponibilidad de las fuentes en GitHub no otorga permisos adicionales de uso comercial o redistribución de recursos ajenos.

## rnai-designer

La capa de integración desarrollada para este proyecto incluye CLI, contratos de datos, automatización headless, análisis de ventanas, reportes, batch, ranking y optimizaciones. Debe distinguirse de los algoritmos originales de si-Fi y de los ejecutables externos utilizados.

**Pendiente de definición por los titulares:** autoría institucional/personal formal y licencia de distribución de `rnai-designer`. No se encontró un archivo de licencia propia ni una declaración de licencia en `pyproject.toml`. Esta edición no asigna MIT, GPL, Creative Commons ni otro permiso al desarrollo propio.

La interacción con componentes licenciados y la distribución conjunta deben revisarse antes de publicar un paquete redistribuible. Este documento no resuelve la compatibilidad jurídica de todas las combinaciones posibles.

## Herramientas, bibliotecas y datos

| Componente | Papel | Tratamiento documental |
|---|---|---|
| Bowtie | Alineamientos de siRNAs | Software externo; versión observada y procedencia a conservar |
| RNAplfold / ViennaRNA | Accesibilidad estructural | Software externo; versión exacta pendiente en esta instalación |
| Python, Conda, Qt/PyQt y dependencias científicas | Entornos de ejecución | Cada distribución conserva sus propias condiciones |
| Matplotlib y dependencias Python 3 | Figuras y soporte del paquete | Instaladas según manifiestos; no son desarrollo propio |
| Índices y FASTA biológicos | Secuencias de referencia | Origen, versión, condiciones y atribuciones por documentar para cada dataset |
| Recursos gráficos/DLL de la copia original | Soporte legado | No asumir que la licencia de una parte cubre todos los recursos |

No se verificaron exhaustivamente las licencias de cada binario, DLL y dataset fuera del repositorio. No se autoriza su redistribución por estar presentes en una instalación local.

## Antes de una distribución pública

Registrar la fuente exacta de si-Fi y cambios locales; conservar créditos y avisos originales; confirmar la licencia del código propio; inventariar las licencias/versiones de binarios y bibliotecas; y documentar origen/condiciones de los índices. Si un dato falta, mantenerlo como pendiente en lugar de inventar una atribución o permiso.

Las rutas locales y los resultados de investigación no son material que deba publicarse automáticamente. `.gitignore` excluye `output/`, pero eso no sustituye una revisión del contenido antes de un commit o publicación.
