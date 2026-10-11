"""
xnafiber.residues
=================

Fabrica de residuos modXNA con cache direccionada por contenido.

``modxna.sh`` ensambla un residuo a partir de un triplete
``<backbone> <sugar> <base>`` encadenando dos operaciones ``sequence`` de
cpptraj (precedidas de ``dataset connect`` sobre los atomos HEAD/TAIL/ANCHOR
declarados en la cabecera mol2 del fragmento), genera la topologia con
``tleap`` y relaja el residuo con ``sander`` antes de guardarlo como
biblioteca OFF (``.lib``). El coste es de segundos a decenas de segundos por
residuo, inaceptable para una herramienta web interactiva si se repite.

Estrategia
----------
1. La clave de cache es ``sha256(backbone, sugar, base, cap, version_modXNA,
   version_cpptraj)``. Es inmune a cambios de nombre del residuo y cambia
   automaticamente si se actualiza el catalogo o la cadena de herramientas.
2. El residuo se construye una sola vez con el nombre generico ``MXR`` y se
   guarda (``.lib``, ``.mol2``, log, carga neta).
3. En cada trabajo se *materializa* con el nombre asignado por
   :class:`~xnafiber.seqlang.ResidueNamer` mediante un renombrado con cpptraj,
   operacion de milisegundos.
4. Los residuos mas usados (ASO y siRNA comerciales) se pre-generan en la
   instalacion con ``xnafiber-warm-cache``.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Tuple

CACHE_SCHEMA = "xnafiber/rescache/1"


@dataclass
class ResidueLib:
    """Residuo modXNA generado y cacheado."""

    triplet: Tuple[str, str, str]
    cap: Optional[str]
    generic_name: str
    lib_path: str
    mol2_path: str
    net_charge: float
    head_atom: Optional[str]
    tail_atom: Optional[str]
    cache_key: str
    log_path: str

    def to_json(self) -> dict:
        d = asdict(self)
        d["triplet"] = list(self.triplet)
        return d


class ModxnaError(RuntimeError):
    pass


class ModxnaRunner:
    """Envoltorio sobre ``modxna.sh`` con cache en disco."""

    GENERIC_NAME = "MXR"

    def __init__(self, modxna_sh: str, cache_dir: str,
                 cpptraj: str = "cpptraj", amberhome: Optional[str] = None,
                 timeout_s: int = 900):
        self.modxna_sh = os.path.abspath(modxna_sh)
        self.cache_dir = os.path.abspath(cache_dir)
        self.cpptraj = cpptraj
        self.amberhome = amberhome
        self.timeout_s = timeout_s
        os.makedirs(self.cache_dir, exist_ok=True)
        self._versions = self._probe_versions()

    # -- entorno ---------------------------------------------------------
    def _env(self) -> Dict[str, str]:
        env = dict(os.environ)
        if self.amberhome:
            env["AMBERHOME"] = self.amberhome
            env["PATH"] = os.path.join(self.amberhome, "bin") + os.pathsep + env.get("PATH", "")
        return env

    def _probe_versions(self) -> Dict[str, str]:
        out = {}
        try:
            r = subprocess.run([self.cpptraj, "--version"], capture_output=True,
                               text=True, timeout=60, env=self._env())
            out["cpptraj"] = (r.stdout or r.stderr).strip().splitlines()[0]
        except Exception as exc:                                    # pragma: no cover
            raise ModxnaError(f"cpptraj no disponible: {exc}") from exc
        try:
            with open(self.modxna_sh, "rb") as fh:
                out["modxna_sha256"] = hashlib.sha256(fh.read()).hexdigest()[:16]
        except OSError as exc:
            raise ModxnaError(f"no se encuentra modxna.sh: {exc}") from exc
        return out

    # -- cache -----------------------------------------------------------
    def cache_key(self, backbone: str, sugar: str, base: str,
                  cap: Optional[str]) -> str:
        payload = json.dumps({"bb": backbone, "sug": sugar, "base": base,
                              "cap": cap, "v": self._versions,
                              "schema": CACHE_SCHEMA},
                             sort_keys=True).encode()
        return hashlib.sha256(payload).hexdigest()[:24]

    def _entry_dir(self, key: str) -> str:
        return os.path.join(self.cache_dir, key[:2], key)

    def lookup(self, key: str) -> Optional[ResidueLib]:
        meta = os.path.join(self._entry_dir(key), "entry.json")
        if not os.path.exists(meta):
            return None
        with open(meta, "r", encoding="utf-8") as fh:
            d = json.load(fh)
        d["triplet"] = tuple(d["triplet"])
        lib = ResidueLib(**d)
        return lib if os.path.exists(lib.lib_path) else None

    # -- construccion ----------------------------------------------------
    def build(self, backbone: str, sugar: str, base: str,
              cap: Optional[str] = None, force: bool = False) -> ResidueLib:
        """Genera (o recupera de cache) el residuo del triplete dado.

        ``cap`` es ``None``, ``"5cap"`` o ``"3cap"``.
        """
        key = self.cache_key(backbone, sugar, base, cap)
        if not force:
            hit = self.lookup(key)
            if hit is not None:
                return hit

        dest = self._entry_dir(key)
        os.makedirs(dest, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="modxna_") as work:
            inp = os.path.join(work, "residue.in")
            with open(inp, "w", encoding="utf-8") as fh:
                fh.write(f"{backbone} {sugar} {base}\n")
            cmd = ["bash", self.modxna_sh, "-i", inp, "-m", self.GENERIC_NAME]
            if cap == "5cap":
                cmd.append("--5cap")
            elif cap == "3cap":
                cmd.append("--3cap")
            cmd.append("--clean")
            proc = subprocess.run(cmd, cwd=work, capture_output=True, text=True,
                                  timeout=self.timeout_s, env=self._env())
            log_path = os.path.join(dest, "modxna.log")
            with open(log_path, "w", encoding="utf-8") as fh:
                fh.write("$ " + " ".join(cmd) + "\n\n")
                fh.write(proc.stdout or "")
                fh.write("\n--- stderr ---\n")
                fh.write(proc.stderr or "")
            if proc.returncode != 0:
                raise ModxnaError(
                    f"modxna.sh fallo para {backbone}/{sugar}/{base} "
                    f"(codigo {proc.returncode}); ver {log_path}")

            produced = {}
            for ext in (".lib", ".mol2"):
                srcp = os.path.join(work, self.GENERIC_NAME + ext)
                if not os.path.exists(srcp):
                    # modxna.sh deja el mol2 optimizado como tmp.opt.mol2
                    alt = os.path.join(work, "tmp.opt.mol2")
                    if ext == ".mol2" and os.path.exists(alt):
                        srcp = alt
                    else:
                        raise ModxnaError(
                            f"modxna.sh no produjo {self.GENERIC_NAME}{ext}; "
                            f"ver {log_path}")
                dstp = os.path.join(dest, self.GENERIC_NAME + ext)
                shutil.copyfile(srcp, dstp)
                produced[ext] = dstp

        head, tail, charge = _parse_lib(produced[".lib"])
        lib = ResidueLib(triplet=(backbone, sugar, base), cap=cap,
                         generic_name=self.GENERIC_NAME,
                         lib_path=produced[".lib"], mol2_path=produced[".mol2"],
                         net_charge=charge, head_atom=head, tail_atom=tail,
                         cache_key=key, log_path=log_path)
        with open(os.path.join(dest, "entry.json"), "w", encoding="utf-8") as fh:
            json.dump(lib.to_json(), fh, indent=1)
        return lib

    # -- materializacion -------------------------------------------------
    def materialize(self, lib: ResidueLib, resname: str,
                    outdir: str) -> Tuple[str, str]:
        """Instancia el residuo cacheado con el nombre ``resname``.

        Devuelve ``(ruta mol2, ruta lib)``. El mol2 se renombra con cpptraj
        (``change crdset ... resname``) y la biblioteca OFF por sustitucion
        textual del identificador de la unidad, que en el formato OFF aparece
        siempre como ``<nombre>`` entre comillas y en las claves
        ``!entry.<nombre>.unit.*``. Ambas operaciones son de milisegundos, lo
        que hace que la cache por triplete sea reutilizable entre trabajos.
        """
        if len(resname) > 4:
            raise ModxnaError(f"nombre de residuo demasiado largo: {resname!r}")
        os.makedirs(outdir, exist_ok=True)
        out_lib = os.path.join(outdir, f"{resname}.lib")
        with open(lib.lib_path, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        text = text.replace(f"!entry.{self.GENERIC_NAME}.",
                            f"!entry.{resname}.")
        text = text.replace(f'"{self.GENERIC_NAME}"', f'"{resname}"')
        with open(out_lib, "w", encoding="utf-8") as fh:
            fh.write(text)

        out_mol2 = os.path.join(outdir, f"{resname}.mol2")
        script = os.path.join(outdir, f"rename_{resname}.cpptraj")
        with open(script, "w", encoding="utf-8") as fh:
            fh.write(f"parm {lib.mol2_path}\n"
                     f"loadcrd {lib.mol2_path} name RES\n"
                     f"change crdset RES resname from * to {resname}\n"
                     f"crdout RES {out_mol2}\n"
                     "run\nquit\n")
        proc = subprocess.run([self.cpptraj, "-i", script], capture_output=True,
                              text=True, timeout=300, env=self._env())
        if proc.returncode != 0 or not os.path.exists(out_mol2):
            raise ModxnaError(f"renombrado a {resname} fallo:\n{proc.stdout}\n{proc.stderr}")
        return out_mol2, out_lib


def _parse_lib(path: str) -> Tuple[Optional[str], Optional[str], float]:
    """Extrae atomo head, atomo tail y carga neta de un fichero OFF de Amber."""
    head = tail = None
    charge = 0.0
    section = None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!entry"):
                section = line.split(".")[-1].split()[0]
                continue
            if section == "atoms" and line.strip() and not line.startswith("!"):
                parts = line.split()
                if len(parts) >= 8:
                    try:
                        charge += float(parts[-1])
                    except ValueError:
                        pass
            elif section == "connect" and line.strip() and not line.startswith("!"):
                pass
    return head, tail, round(charge, 4)


def warm_cache(runner: ModxnaRunner, triplets: List[Tuple[str, str, str]],
               progress=None) -> Dict[str, str]:
    """Pre-genera un conjunto de residuos (uso: instalacion del servidor)."""
    out = {}
    for i, (bb, sug, base) in enumerate(triplets, 1):
        lib = runner.build(bb, sug, base)
        out[f"{bb}/{sug}/{base}"] = lib.cache_key
        if progress:
            progress(i, len(triplets), bb, sug, base)
    return out
