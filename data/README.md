# Carpeta de datos (vacía a propósito)

Los datos de ADNI **no se distribuyen** con este repositorio: el Data Use Agreement
de ADNI lo prohíbe. Deben obtenerse en https://ida.loni.usc.edu previa solicitud.

Colocar aquí los siguientes ficheros y ejecutar `src/01_merge.py`, que generará
`data/adni_merged.csv` (630 × 22):

| Fichero | Ruta en la web de ADNI |
|---|---|
| `ADNIMERGE2_*.tar.gz` | Study Data → Study Info → Data & Database |
| `UWNPSYCHSUM_*.csv` | Assessments → Neuropsychological |
| `UCSFFSX51_11_08_19_*.csv` | Imaging → MR Image Analysis |
| `UCBERKELEY_AMY_6MM_*.csv` | Imaging → PET Image Analysis |
| `UCBERKELEYFDG_8mm_02_17_23_*.csv` | Imaging → PET Image Analysis |
| `DATADIC_*.csv` | Study Data → Study Info → Data & Database |

Los sufijos de fecha varían según la descarga. Anotar la fecha exacta: se declara
en el artículo.
