"""Gráficos descriptivos exportables, sin interpolar observaciones faltantes."""
from pathlib import Path
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg


def write_plot(path: Path, points: list[dict], gene_id: str, length: int,
               simulated: bool = False, evaluated: bool = True) -> None:
    figure = Figure(figsize=(10, 5), layout="constrained")
    FigureCanvasAgg(figure)
    access_ax, label_ax = figure.subplots(2, 1, sharex=True)
    colors = {True: "#16853b", False: "#c83232", None: "#777777"}
    for label, color in colors.items():
        group = [point for point in points if point["is_efficient"] is label]
        known = [point for point in group if point["accessibility_value"] is not None]
        access_ax.scatter([p["start_1based"] for p in known],
                          [p["accessibility_value"] for p in known], color=color)
        label_ax.scatter([p["start_1based"] for p in group],
                         [{True: 1, False: 0, None: -1}[label]] * len(group), color=color)
    missing = [p for p in points if p["accessibility_value"] is None]
    access_ax.scatter([p["start_1based"] for p in missing], [-0.08] * len(missing),
                      marker="x", color="#777777")
    access_ax.set_ylim(-0.15, 1.05)
    access_ax.set_ylabel("accessibility_value")
    access_ax.set_title("Accesibilidad importada; cruces bajo cero: dato faltante/conflictivo")
    label_ax.set_yticks([-1, 0, 1], ["null/conflicto", "false", "true"])
    label_ax.set_ylim(-1.5, 1.5)
    label_ax.set_ylabel("is_efficient (si-Fi)")
    label_ax.set_xlabel("Inicio del siRNA sobre la referencia (1-based)")
    label_ax.set_xlim(0.5, max(1, length) + 0.5)
    for axis in (access_ax, label_ax):
        axis.grid(alpha=0.2)
    title = f"{gene_id} — métricas de si-Fi"
    if not evaluated:
        title += "\nRNAi Designer: NO EVALUADO"
    if simulated:
        title += "\nDATOS FICTICIOS — sin ejecución real"
    figure.suptitle(title)
    figure.supxlabel("Predicción computacional: no demuestra eficacia ni seguridad", fontsize=9)
    if not points:
        access_ax.text(0.5, 0.5, "Sin observaciones: no permite inferir ausencia de riesgo",
                       transform=access_ax.transAxes, ha="center")
    figure.savefig(path, format="svg", metadata={"Date": None})
    figure.savefig(path.with_suffix(".png"), dpi=160)
    figure.clear()


def write_offtarget_plot(path: Path, evaluations: list[dict], gene_id: str, length: int,
                         simulated: bool = False) -> None:
    figure = Figure(figsize=(12, 5), layout="constrained")
    FigureCanvasAgg(figure)
    axis = figure.subplots()
    labels = []
    for y, entry in enumerate(evaluations):
        criterion = entry["criterion"]
        labels.append(f'{entry["organism"]}\n{criterion["label"]}')
        for match_type, color, marker in ((True, "#bd2727", "o"), (False, "#3974ad", "."), (None, "#777777", "x")):
            group = [a for a in criterion["alignments"] if a["contains_exact_21nt"] is match_type]
            axis.scatter([a["query_start_1based"] for a in group], [y] * len(group), color=color, marker=marker)
        if not criterion["alignments"]:
            axis.text(0.02, y, "Sin filas de alineamiento; interpretar según estado", va="center",
                      transform=axis.get_yaxis_transform(), fontsize=9)
    axis.set_yticks(range(len(evaluations)), labels)
    axis.set_ylim(-0.7, max(0, len(evaluations)-1) + 0.7)
    axis.invert_yaxis()
    axis.set_xlim(0.5, max(1, length) + 0.5)
    axis.set_xlabel("Inicio en query (1-based). Rojo: contiene 21nt exactos; azul: no; gris: sin verificar")
    axis.grid(alpha=0.2)
    title = f"{gene_id} — coincidencias importadas por especie"
    if simulated:
        title += "\nQUERY/DATOS DE DEMOSTRACIÓN — sin validación experimental"
    axis.set_title(title)
    figure.supxlabel("El criterio se limita a las bases consultadas; no garantiza seguridad ambiental", fontsize=9)
    figure.savefig(path, format="svg", metadata={"Date": None})
    figure.savefig(path.with_suffix(".png"), dpi=160)
    figure.clear()
