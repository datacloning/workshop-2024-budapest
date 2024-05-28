# Szükséges szoftver telepítése

## JAGS

Követsd az instrukciókat a JAGS honlapon <https://mcmc-jags.sourceforge.io/>.

```bash
# Ubuntu Linux
sudo apt update -qq && sudo apt install --yes --no-install-recommends jags

# Mac
brew install jags
```

Windows: töltsd le [innen](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/).

## R csomagok

R-ben:

```R
# CRAN packages
install.packages(c("rjags", "shiny", "dclone", "R2WinBUGS"))
```

## Teszt

Ha a következő példa le tud futni, akkor minden sikeresen telepítve:

```R
example("jags.fit", package = "dclone", run.dontrun = TRUE)
```
