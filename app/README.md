# Shiny applikáció

Az applikáció segítségével különféle statisztikai eloszlásokat próbálhatunk ki, vagy megismerkedhetünk a Bernoulli modellel mind frekventista, mind Bayes-i szempontból. Megtapasztalhatjuk a prior eloszlás hatását, és megláthatjuk hogyan működik az adat klónozás.

A futtatás ebből a könyvtárból R-ben a következő parancs segítségével: `shiny::runApp()`.

A már telepített verzió itt érhető el: <https://psolymos.shinyapps.io/dcapps/>.

Shinylive webr segítségével (<https://posit-dev.github.io/r-shinylive/>):

```R
shinylive::export("app", "docs")
httpuv::runStaticServer("docs/", port=8080)
```

A Shinylive verzió itt érhető el: <https://datacloning.org/workshop-2024-budapest/>.
