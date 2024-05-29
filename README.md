# Adat Klónozás Kurzus - Data Cloning Workshop / Budapest 2014

> A kurzus nyelve **angol**!

## Részvétel (MÉG ALAKULÓBAN)

A kurzus 8 óra hosszú, 2024 augusztus 28-án, de. 8-tól du. 5-ig, az [Állatorvostudományi Egyetem](https://maps.app.goo.gl/AaZDqDmGhdfs3n4r9) [Biostatisztika Tanszékének](https://univet.hu/hu/egyetem/szervezeti-egysegek/biomatematikai-es-szamitastechnikai-tanszek/) számítógépes laborjában (1078 Budapest, István utca 2., N épület, 3. emelet, N3 terem).

Jelentkezés az [online űrlap kitöltésével](https://forms.reform.app/analythium/data-cloning-workshop-2024/npne8q) lehetséges.

## Oktatók

- [Subhash Lele](https://scholar.google.ca/citations?hl=en&user=1CNJm5UAAAAJ­), Professor Emeritus ([Dept. of Mathematical and Statistical Sciences, University of Alberta, Edmonton](https://sites.ualberta.ca/~slele/))
- [Peter Solymos](https://peter.solymos.org/­), szenior adattudós ([E Source](https://esource.com)), címzetes egyetemi tanár ([Dept. of Biological Sciences, University of Alberta, Edmonton](https://www.ualberta.ca/biological-sciences/faculty-and-staff/lecturers-adjunct/index.html))

## Összefoglaló

A kevert, vagy hierarchikus modellek nagyon hasznosak sokféle alkalmazott területen. A kurzus célja, hogy különféle gyakorlati problémák bemutatása révén bevezetést nyújtson ezen modellek felépítésébe, elméletébe, és implementációjába.

A kurzus során előadások és számítógépes gyakorlatok révén fogjuk elmagyarázni a legyakoribb modellek felépítését és működését. Ezt aztán a résztvevők saját igényeiknek megfelelően tovább fejleszthetik.

A résztvevőktől azt várjuk, hogy ismerik az R programozási nyelv alapjait és jártasak a regressziószámítás területén. A hierarchikus modellek vagy a Bayes-i statisztika ismerete nem követelmény, de hasznos lehet.

A kurzus végén a résztvevők megfelelő alapokkal fognak rendelkezni ahhoz, hogy kritikusan tudjanak gondolkozni a frekventista és Bayes-i modellekről, képesek legyenek (majdnem) minden hierarchikus modell használatára.


## Előkészületek

Követsd az [instrukciókat](setup.md) és hozd magaddal a laptopodat.

Próbáld ki a [Shiny applikációt](./app/).

## A kurzus alatt

Töltsd le a [bezippelt](https://github.com/datacloning/workshop-2024-budapest/archive/refs/heads/main.zip) verzióját ennek a repónak, vagy használd a `git clone https://github.com/datacloning/workshop-2023-edmonton.git` utasítást ha jártas vagy a Git-ben.

## Jegyzetek (FRISSÍTÉS ALATT)

A jegyzetek markdown formátumban a GitHub-on érhetőek el, ajánljuk a _világos_ mód használatát, hogy a matematikai képletek jól látszádjanak.

| Témakör    | Linkek |
| -------- | ------- |
| Előkészületek  | [Előkészületek](setup.md)  |
| Part 1. Bevezetés  | [Jegyzetek](./01-intro/), [App](./app/)  |
| Part 2. Kevert modellek  | [Jegyzetek](./02-mixed-models/)  |
| Part 3. Idősorok  | [Jegyzetek](./03-time-series/)  |
| Part 4. Gyakorlati magfontolások  | [Jegyzetek](./04-other/)  |
| Irodalmak  | [PDF fájlok](./docs/)  |

## Licensz

A kurzus licensze az
[Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)](https://creativecommons.org/licenses/by-sa/4.0/).
Publikációkban kérjük az alábbi cikkeket hivatkozni az elmélet és a szofveres implementációra vonatkozóan:

- Lele, S. R., and Solymos, P., 2023. Data Cloning Workshop - Hierarchical Models Made Easy. April 13, 2023. URL <https://github.com/datacloning/workshop-2023-edmonton>
- Lele, S.R., B. Dennis and F. Lutscher, 2007. Data cloning: easy maximum likelihood estimation for complex ecological models using Bayesian Markov chain Monte Carlo methods. Ecology Letters 10, 551-563. [DOI 10.1111/j.1461-0248.2007.01047.x­](https://doi.org/10.1111/j.1461-0248.2007.01047.x)
- Lele, S. R., Nadeem, K., and Schmuland, B., 2010. Estimability and likelihood inference for generalized linear mixed models using data cloning. Journal of the American Statistical Association 105, 1617-1625. [DOI 10.1198/jasa.2010.tm09757­](https://doi.org/10.1198/jasa.2010.tm09757)
- Solymos, P., 2010. dclone: Data Cloning in R. The R Journal 2(2), 29-37. URL <https://journal.r-project.org/archive/2010/RJ-2010-011/RJ-2010-011.pdf>

A felhasznált szofverek (JAGS, rjags, R2WinBUGS, dclone) licensze [GPL-2](https://cran.r-project.org/web/licenses/GPL-2).
