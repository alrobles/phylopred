# Prepare an incidence matrix from an interaction table

Builds an incidence (presence/absence) matrix from a two-column
data.frame of species interactions.

## Usage

``` r
prepare_incidence_matrix(db)
```

## Arguments

- db:

  A data.frame with exactly two columns (e.g., parasite and host).

## Value

An incidence matrix (rows = first column unique values, cols = second
column unique values) filled with 0/1.

## Examples

``` r
prepare_incidence_matrix(beetleTreeInteractions)
#>                           Lysiloma Lindera Castanea Trema Cedrela Bursera Ulmus
#> Ambrosiodmus obliquus            1       1        1     1       1       1     0
#> Ambrosiodmus rubricollis         0       0        0     0       0       0     1
#> Ambrosiodmus tachygraphus        0       0        1     0       0       0     0
#> Xyleborus atratus                0       0        0     0       0       0     1
#> Xyleborus obesus                 0       0        0     0       0       0     0
#> mutillatus                       0       1        0     0       0       0     0
#> Coptoborus pseudotenuis          0       0        0     0       1       1     0
#> Euwallacea fornicatus            0       0        0     0       1       0     0
#> Euwallacea validus               0       0        1     0       0       0     1
#> Sampsonius dampfi                0       0        0     0       0       0     0
#> Xyleborinus saxeseni             0       0        0     0       1       0     0
#> Xyleborus affinis                0       0        1     0       1       1     0
#> Xyleborus bispinatus             0       0        0     0       0       0     0
#> Xyleborus ferrugineus            1       0        0     0       1       1     0
#> Xyleborus glabratus              0       1        0     0       0       0     0
#> Xyleborus posticus               0       0        0     0       0       0     0
#> Xyleborus seriatus               0       0        0     0       0       0     0
#> Xyleborus spathipennis           0       0        0     0       0       1     0
#> Xyleborus xylographus            0       0        1     0       0       0     0
#> Xylosandrus compactus            0       0        0     0       0       0     0
#> Xylosandrus crassiusculus        0       0        0     0       0       0     1
#> Xyleborus germanus               0       1        1     0       0       0     1
#> Xylosandrus morigerus            0       0        0     1       0       1     0
#>                           Gardenia Carya Cupressus Terminalia Juglans Prunus
#> Ambrosiodmus obliquus            0     0         0          0       0      0
#> Ambrosiodmus rubricollis         1     1         1          1       1      1
#> Ambrosiodmus tachygraphus        0     1         0          0       0      1
#> Xyleborus atratus                0     0         0          0       0      0
#> Xyleborus obesus                 0     0         0          0       0      0
#> mutillatus                       0     0         0          0       0      0
#> Coptoborus pseudotenuis          0     0         0          0       0      0
#> Euwallacea fornicatus            0     0         0          1       0      0
#> Euwallacea validus               0     0         0          0       0      1
#> Sampsonius dampfi                0     0         0          0       0      0
#> Xyleborinus saxeseni             0     1         0          0       0      1
#> Xyleborus affinis                0     1         0          1       0      1
#> Xyleborus bispinatus             0     0         0          0       0      0
#> Xyleborus ferrugineus            0     0         0          1       0      0
#> Xyleborus glabratus              0     0         0          0       0      0
#> Xyleborus posticus               0     0         1          0       0      0
#> Xyleborus seriatus               0     0         0          0       0      0
#> Xyleborus spathipennis           0     0         0          0       0      0
#> Xyleborus xylographus            0     0         0          0       0      0
#> Xylosandrus compactus            0     1         0          0       0      0
#> Xylosandrus crassiusculus        0     1         0          0       0      1
#> Xyleborus germanus               0     1         0          0       1      0
#> Xylosandrus morigerus            0     0         0          1       0      0
#>                           Liriodendron Alnus Betula Pinus Fagus Populus Acer
#> Ambrosiodmus obliquus                0     0      0     0     0       0    0
#> Ambrosiodmus rubricollis             0     0      0     0     0       0    0
#> Ambrosiodmus tachygraphus            1     0      0     0     0       0    0
#> Xyleborus atratus                    0     1      1     1     0       0    0
#> Xyleborus obesus                     1     0      0     0     1       1    1
#> mutillatus                           0     0      0     0     1       0    0
#> Coptoborus pseudotenuis              0     0      0     0     0       0    0
#> Euwallacea fornicatus                0     0      0     0     0       0    0
#> Euwallacea validus                   0     0      0     1     1       1    0
#> Sampsonius dampfi                    0     0      0     0     0       0    0
#> Xyleborinus saxeseni                 1     0      1     1     0       1    1
#> Xyleborus affinis                    0     0      1     1     0       1    1
#> Xyleborus bispinatus                 0     0      0     0     0       0    0
#> Xyleborus ferrugineus                0     0      0     1     0       0    1
#> Xyleborus glabratus                  0     0      0     0     0       0    0
#> Xyleborus posticus                   0     0      0     0     0       0    0
#> Xyleborus seriatus                   0     0      0     0     1       0    1
#> Xyleborus spathipennis               0     0      0     0     0       0    0
#> Xyleborus xylographus                0     0      0     0     0       0    0
#> Xylosandrus compactus                0     0      0     0     0       0    1
#> Xylosandrus crassiusculus            1     0      0     1     0       1    0
#> Xyleborus germanus                   1     0      1     1     1       0    1
#> Xylosandrus morigerus                0     0      0     0     0       0    0
#>                           Ormosia Cinnamomum Cercis Osmanthus Cryptomeria
#> Ambrosiodmus obliquus           0          0      0         0           0
#> Ambrosiodmus rubricollis        0          0      0         0           0
#> Ambrosiodmus tachygraphus       0          0      0         0           0
#> Xyleborus atratus               0          0      0         0           0
#> Xyleborus obesus                0          0      0         0           0
#> mutillatus                      1          1      1         1           1
#> Coptoborus pseudotenuis         0          0      0         0           0
#> Euwallacea fornicatus           0          0      0         0           0
#> Euwallacea validus              0          0      0         0           1
#> Sampsonius dampfi               0          0      0         0           0
#> Xyleborinus saxeseni            0          0      0         0           0
#> Xyleborus affinis               0          0      0         0           0
#> Xyleborus bispinatus            0          0      0         0           0
#> Xyleborus ferrugineus           0          0      0         0           0
#> Xyleborus glabratus             0          1      0         0           0
#> Xyleborus posticus              0          0      0         0           0
#> Xyleborus seriatus              0          0      0         0           1
#> Xyleborus spathipennis          0          0      0         0           0
#> Xyleborus xylographus           0          0      0         0           0
#> Xylosandrus compactus           0          1      0         0           0
#> Xylosandrus crassiusculus       0          0      1         0           0
#> Xyleborus germanus              0          0      0         0           0
#> Xylosandrus morigerus           0          0      0         0           0
#>                           Carpinus Swietenia Toxicodendron Machilus Artocarpus
#> Ambrosiodmus obliquus            0         0             0        0          0
#> Ambrosiodmus rubricollis         0         0             0        0          0
#> Ambrosiodmus tachygraphus        0         0             0        0          0
#> Xyleborus atratus                0         0             0        0          0
#> Xyleborus obesus                 0         0             0        0          0
#> mutillatus                       1         1             1        1          0
#> Coptoborus pseudotenuis          0         0             0        0          1
#> Euwallacea fornicatus            0         0             0        0          0
#> Euwallacea validus               1         0             0        0          0
#> Sampsonius dampfi                0         0             0        0          0
#> Xyleborinus saxeseni             0         0             0        0          0
#> Xyleborus affinis                0         1             1        0          0
#> Xyleborus bispinatus             0         1             0        0          0
#> Xyleborus ferrugineus            0         1             0        0          1
#> Xyleborus glabratus              0         0             0        1          0
#> Xyleborus posticus               0         0             0        0          0
#> Xyleborus seriatus               1         0             0        0          0
#> Xyleborus spathipennis           0         0             0        0          0
#> Xyleborus xylographus            0         0             0        0          0
#> Xylosandrus compactus            0         1             0        0          0
#> Xylosandrus crassiusculus        0         1             0        1          0
#> Xyleborus germanus               1         0             1        0          0
#> Xylosandrus morigerus            0         1             0        0          0
#>                           Heliocarpus Hevea Theobroma Dialium Vachellia Delonix
#> Ambrosiodmus obliquus               0     0         0       0         0       0
#> Ambrosiodmus rubricollis            0     0         0       0         0       0
#> Ambrosiodmus tachygraphus           0     0         0       0         0       0
#> Xyleborus atratus                   0     0         0       0         0       0
#> Xyleborus obesus                    0     0         0       0         0       0
#> mutillatus                          0     0         0       0         0       0
#> Coptoborus pseudotenuis             1     1         1       1         1       1
#> Euwallacea fornicatus               0     1         1       0         0       1
#> Euwallacea validus                  0     0         0       0         0       0
#> Sampsonius dampfi                   0     0         1       0         0       0
#> Xyleborinus saxeseni                0     0         0       0         0       0
#> Xyleborus affinis                   0     0         1       1         0       0
#> Xyleborus bispinatus                0     0         0       0         0       0
#> Xyleborus ferrugineus               1     0         1       1         0       1
#> Xyleborus glabratus                 0     0         0       0         0       0
#> Xyleborus posticus                  1     0         1       0         0       0
#> Xyleborus seriatus                  0     0         0       0         0       0
#> Xyleborus spathipennis              0     0         0       0         0       0
#> Xyleborus xylographus               0     0         0       0         0       0
#> Xylosandrus compactus               0     0         0       0         0       0
#> Xylosandrus crassiusculus           0     1         1       0         0       0
#> Xyleborus germanus                  0     0         0       0         0       0
#> Xylosandrus morigerus               0     0         1       0         0       0
#>                           Gliricidia Citrus Persea Gmelina Pterocymbium
#> Ambrosiodmus obliquus              0      0      0       0            0
#> Ambrosiodmus rubricollis           0      0      0       0            0
#> Ambrosiodmus tachygraphus          0      0      0       0            0
#> Xyleborus atratus                  0      0      0       0            0
#> Xyleborus obesus                   0      0      0       0            0
#> mutillatus                         0      0      0       0            0
#> Coptoborus pseudotenuis            1      1      0       0            0
#> Euwallacea fornicatus              0      0      1       1            1
#> Euwallacea validus                 0      0      0       0            0
#> Sampsonius dampfi                  0      0      0       0            0
#> Xyleborinus saxeseni               0      0      0       0            0
#> Xyleborus affinis                  0      0      1       0            0
#> Xyleborus bispinatus               0      0      1       0            0
#> Xyleborus ferrugineus              1      0      1       0            0
#> Xyleborus glabratus                0      0      1       0            0
#> Xyleborus posticus                 0      0      0       0            0
#> Xyleborus seriatus                 0      0      0       0            0
#> Xyleborus spathipennis             0      0      0       0            0
#> Xyleborus xylographus              0      0      0       0            0
#> Xylosandrus compactus              0      0      1       0            0
#> Xylosandrus crassiusculus          1      0      0       0            0
#> Xyleborus germanus                 0      0      0       0            0
#> Xylosandrus morigerus              0      0      1       0            0
#>                           Tephrosia Bombax Ricinus Casuarina Milicia Falcataria
#> Ambrosiodmus obliquus             0      0       0         0       0          0
#> Ambrosiodmus rubricollis          0      0       0         0       0          0
#> Ambrosiodmus tachygraphus         0      0       0         0       0          0
#> Xyleborus atratus                 0      0       0         0       0          0
#> Xyleborus obesus                  0      0       0         0       0          0
#> mutillatus                        0      0       0         0       0          0
#> Coptoborus pseudotenuis           0      0       0         0       0          0
#> Euwallacea fornicatus             1      1       1         1       1          1
#> Euwallacea validus                0      0       0         0       0          0
#> Sampsonius dampfi                 0      0       0         0       0          0
#> Xyleborinus saxeseni              0      0       0         1       0          0
#> Xyleborus affinis                 0      0       0         0       0          0
#> Xyleborus bispinatus              0      0       0         0       0          0
#> Xyleborus ferrugineus             0      0       0         0       0          0
#> Xyleborus glabratus               0      0       0         0       0          0
#> Xyleborus posticus                0      0       0         0       0          0
#> Xyleborus seriatus                0      0       0         0       0          0
#> Xyleborus spathipennis            0      0       0         0       0          0
#> Xyleborus xylographus             0      0       0         0       0          0
#> Xylosandrus compactus             0      0       0         0       0          0
#> Xylosandrus crassiusculus         0      0       0         0       0          0
#> Xyleborus germanus                0      0       0         0       0          0
#> Xylosandrus morigerus             1      0       0         0       0          1
#>                           Myristica Mangifera Clerodendrum Ochroma Erythrina
#> Ambrosiodmus obliquus             0         0            0       0         0
#> Ambrosiodmus rubricollis          0         0            0       0         0
#> Ambrosiodmus tachygraphus         0         0            0       0         0
#> Xyleborus atratus                 0         0            0       0         0
#> Xyleborus obesus                  0         0            0       0         0
#> mutillatus                        0         0            0       0         0
#> Coptoborus pseudotenuis           0         0            0       0         0
#> Euwallacea fornicatus             1         1            1       1         1
#> Euwallacea validus                0         0            0       0         0
#> Sampsonius dampfi                 0         0            0       0         0
#> Xyleborinus saxeseni              0         0            0       0         0
#> Xyleborus affinis                 0         1            0       0         1
#> Xyleborus bispinatus              0         0            0       0         0
#> Xyleborus ferrugineus             0         1            0       0         1
#> Xyleborus glabratus               0         0            0       0         0
#> Xyleborus posticus                0         0            0       0         1
#> Xyleborus seriatus                0         0            0       0         0
#> Xyleborus spathipennis            0         0            0       0         0
#> Xyleborus xylographus             0         0            0       0         0
#> Xylosandrus compactus             0         0            0       0         0
#> Xylosandrus crassiusculus         0         0            0       0         0
#> Xyleborus germanus                0         0            0       0         0
#> Xylosandrus morigerus             0         1            0       1         0
#>                           Albizia Protium Ixora Robinia Callerya Grevillea
#> Ambrosiodmus obliquus           0       0     0       0        0         0
#> Ambrosiodmus rubricollis        0       0     0       0        0         0
#> Ambrosiodmus tachygraphus       0       0     0       0        0         0
#> Xyleborus atratus               0       0     0       0        0         0
#> Xyleborus obesus                0       0     0       0        0         0
#> mutillatus                      0       0     0       0        0         0
#> Coptoborus pseudotenuis         0       0     0       0        0         0
#> Euwallacea fornicatus           1       1     1       1        1         1
#> Euwallacea validus              0       0     0       0        0         0
#> Sampsonius dampfi               0       0     0       0        0         0
#> Xyleborinus saxeseni            0       0     0       0        0         0
#> Xyleborus affinis               1       0     0       0        0         0
#> Xyleborus bispinatus            0       0     0       0        0         0
#> Xyleborus ferrugineus           0       0     0       0        0         0
#> Xyleborus glabratus             0       0     0       0        0         0
#> Xyleborus posticus              0       0     0       0        0         0
#> Xyleborus seriatus              0       0     0       0        0         0
#> Xyleborus spathipennis          0       0     0       0        0         0
#> Xyleborus xylographus           0       0     0       0        0         0
#> Xylosandrus compactus           0       0     0       0        0         0
#> Xylosandrus crassiusculus       1       0     0       0        0         1
#> Xyleborus germanus              0       0     0       0        0         0
#> Xylosandrus morigerus           0       0     0       0        0         0
#>                           Shorea Polyscias Senna Camellia Crotalaria Lannea
#> Ambrosiodmus obliquus          0         0     0        0          0      0
#> Ambrosiodmus rubricollis       0         0     0        0          0      0
#> Ambrosiodmus tachygraphus      0         0     0        0          0      0
#> Xyleborus atratus              0         0     0        0          0      0
#> Xyleborus obesus               0         0     0        0          0      0
#> mutillatus                     0         0     0        0          0      0
#> Coptoborus pseudotenuis        0         0     0        0          0      0
#> Euwallacea fornicatus          1         1     1        1          1      1
#> Euwallacea validus             0         0     0        0          0      0
#> Sampsonius dampfi              0         0     0        0          0      0
#> Xyleborinus saxeseni           0         0     0        0          0      0
#> Xyleborus affinis              0         0     0        0          0      0
#> Xyleborus bispinatus           0         0     0        0          0      0
#> Xyleborus ferrugineus          0         0     0        0          0      0
#> Xyleborus glabratus            1         0     0        0          0      0
#> Xyleborus posticus             0         0     0        0          0      0
#> Xyleborus seriatus             0         0     0        0          0      0
#> Xyleborus spathipennis         0         0     0        0          0      0
#> Xyleborus xylographus          0         0     0        0          0      0
#> Xylosandrus compactus          0         0     0        1          0      0
#> Xylosandrus crassiusculus      1         0     0        0          0      1
#> Xyleborus germanus             0         0     1        0          0      0
#> Xylosandrus morigerus          0         0     1        1          1      0
#>                           Xylia Ailanthus Phellodendron Tilia Aphananthe Ficus
#> Ambrosiodmus obliquus         0         0             0     0          0     0
#> Ambrosiodmus rubricollis      0         0             0     0          0     0
#> Ambrosiodmus tachygraphus     0         0             0     0          0     0
#> Xyleborus atratus             0         0             0     0          0     0
#> Xyleborus obesus              0         0             0     0          0     0
#> mutillatus                    0         0             0     0          0     0
#> Coptoborus pseudotenuis       0         0             0     0          0     0
#> Euwallacea fornicatus         1         0             0     0          0     0
#> Euwallacea validus            0         1             1     1          1     1
#> Sampsonius dampfi             0         0             0     0          0     0
#> Xyleborinus saxeseni          0         0             0     1          0     0
#> Xyleborus affinis             0         0             0     0          0     0
#> Xyleborus bispinatus          0         0             0     0          0     0
#> Xyleborus ferrugineus         0         0             0     0          0     1
#> Xyleborus glabratus           0         0             0     0          0     0
#> Xyleborus posticus            0         0             0     0          0     1
#> Xyleborus seriatus            0         0             0     1          0     0
#> Xyleborus spathipennis        0         0             0     0          0     0
#> Xyleborus xylographus         0         0             0     0          0     0
#> Xylosandrus compactus         0         0             0     0          0     1
#> Xylosandrus crassiusculus     0         0             0     0          0     1
#> Xyleborus germanus            0         0             0     1          0     0
#> Xylosandrus morigerus         0         0             0     0          0     0
#>                           Abies Quercus Dalbergia Cleyera Mallotus Cunninghamia
#> Ambrosiodmus obliquus         0       0         0       0        0            0
#> Ambrosiodmus rubricollis      0       0         0       0        0            0
#> Ambrosiodmus tachygraphus     0       0         0       0        0            0
#> Xyleborus atratus             0       0         0       0        0            0
#> Xyleborus obesus              0       0         0       0        0            0
#> mutillatus                    0       0         0       0        0            0
#> Coptoborus pseudotenuis       0       0         0       0        0            0
#> Euwallacea fornicatus         0       0         0       0        0            0
#> Euwallacea validus            1       1         1       1        1            1
#> Sampsonius dampfi             0       0         0       0        0            0
#> Xyleborinus saxeseni          0       1         0       0        0            0
#> Xyleborus affinis             0       1         0       0        0            0
#> Xyleborus bispinatus          0       0         0       0        0            0
#> Xyleborus ferrugineus         0       1         0       0        0            0
#> Xyleborus glabratus           0       0         0       0        0            0
#> Xyleborus posticus            0       0         0       0        0            0
#> Xyleborus seriatus            0       0         0       1        1            0
#> Xyleborus spathipennis        0       0         0       0        0            0
#> Xyleborus xylographus         0       1         0       0        0            0
#> Xylosandrus compactus         0       1         0       0        0            0
#> Xylosandrus crassiusculus     0       1         1       0        0            0
#> Xyleborus germanus            1       1         0       1        0            0
#> Xylosandrus morigerus         0       0         1       0        0            0
#>                           Magnolia Chamaecyparis Zelkova Tsuga Celtis Miconia
#> Ambrosiodmus obliquus            0             0       0     0      0       0
#> Ambrosiodmus rubricollis         0             0       0     0      0       0
#> Ambrosiodmus tachygraphus        0             0       0     0      0       0
#> Xyleborus atratus                0             0       0     0      0       0
#> Xyleborus obesus                 0             0       0     0      0       0
#> mutillatus                       0             0       0     0      0       0
#> Coptoborus pseudotenuis          0             0       0     0      0       0
#> Euwallacea fornicatus            0             0       0     0      0       0
#> Euwallacea validus               1             1       1     1      1       0
#> Sampsonius dampfi                0             0       0     0      0       1
#> Xyleborinus saxeseni             0             0       0     0      1       0
#> Xyleborus affinis                0             1       0     0      1       1
#> Xyleborus bispinatus             0             0       0     0      0       0
#> Xyleborus ferrugineus            0             0       0     0      0       1
#> Xyleborus glabratus              0             0       0     0      0       0
#> Xyleborus posticus               0             0       0     0      0       0
#> Xyleborus seriatus               0             0       0     0      0       0
#> Xyleborus spathipennis           0             0       0     0      0       0
#> Xyleborus xylographus            0             0       0     0      0       0
#> Xylosandrus compactus            1             0       0     0      1       0
#> Xylosandrus crassiusculus        1             0       0     0      0       0
#> Xyleborus germanus               0             1       0     1      1       0
#> Xylosandrus morigerus            0             0       0     0      0       1
#>                           Senegalia Qualea Sorbus Malus Pyrus Actinidia
#> Ambrosiodmus obliquus             0      0      0     0     0         0
#> Ambrosiodmus rubricollis          0      0      0     0     0         0
#> Ambrosiodmus tachygraphus         0      0      0     0     0         0
#> Xyleborus atratus                 0      0      0     0     0         0
#> Xyleborus obesus                  0      0      0     0     0         0
#> mutillatus                        0      0      0     0     0         0
#> Coptoborus pseudotenuis           0      0      0     0     0         0
#> Euwallacea fornicatus             0      0      0     0     0         0
#> Euwallacea validus                0      0      0     0     0         0
#> Sampsonius dampfi                 1      1      0     0     0         0
#> Xyleborinus saxeseni              0      0      1     1     1         1
#> Xyleborus affinis                 1      0      0     0     0         0
#> Xyleborus bispinatus              0      0      0     0     0         0
#> Xyleborus ferrugineus             1      0      0     0     0         0
#> Xyleborus glabratus               0      0      0     0     0         0
#> Xyleborus posticus                0      0      0     0     0         0
#> Xyleborus seriatus                0      0      0     0     0         0
#> Xyleborus spathipennis            0      0      0     0     0         0
#> Xyleborus xylographus             0      0      0     0     0         0
#> Xylosandrus compactus             0      0      0     0     0         0
#> Xylosandrus crassiusculus         0      0      1     1     0         0
#> Xyleborus germanus                0      0      0     0     0         0
#> Xylosandrus morigerus             0      0      0     0     0         0
#>                           Notholithocarpus Taxodium Cornus Eucalyptus Arbutus
#> Ambrosiodmus obliquus                    0        0      0          0       0
#> Ambrosiodmus rubricollis                 0        0      0          0       0
#> Ambrosiodmus tachygraphus                0        0      0          0       0
#> Xyleborus atratus                        0        0      0          0       0
#> Xyleborus obesus                         0        0      0          0       0
#> mutillatus                               0        0      0          0       0
#> Coptoborus pseudotenuis                  0        0      0          0       0
#> Euwallacea fornicatus                    0        0      0          0       0
#> Euwallacea validus                       0        0      0          0       0
#> Sampsonius dampfi                        0        0      0          0       0
#> Xyleborinus saxeseni                     1        1      1          1       1
#> Xyleborus affinis                        0        1      0          1       0
#> Xyleborus bispinatus                     0        0      0          0       0
#> Xyleborus ferrugineus                    0        0      0          1       0
#> Xyleborus glabratus                      0        0      0          0       0
#> Xyleborus posticus                       0        0      0          0       0
#> Xyleborus seriatus                       0        0      0          0       0
#> Xyleborus spathipennis                   0        0      0          0       0
#> Xyleborus xylographus                    0        0      1          0       0
#> Xylosandrus compactus                    0        0      1          0       0
#> Xylosandrus crassiusculus                0        0      1          0       0
#> Xyleborus germanus                       0        1      1          0       0
#> Xylosandrus morigerus                    0        0      0          0       0
#>                           Pseudotsuga Plumeria Annona Diospyros Rhus Brosimum
#> Ambrosiodmus obliquus               0        0      0         0    0        0
#> Ambrosiodmus rubricollis            0        0      0         0    0        0
#> Ambrosiodmus tachygraphus           0        0      0         0    0        0
#> Xyleborus atratus                   0        0      0         0    0        0
#> Xyleborus obesus                    0        0      0         0    0        0
#> mutillatus                          0        0      0         0    0        0
#> Coptoborus pseudotenuis             0        0      0         0    0        0
#> Euwallacea fornicatus               0        0      0         0    0        0
#> Euwallacea validus                  0        0      0         0    0        0
#> Sampsonius dampfi                   0        0      0         0    0        0
#> Xyleborinus saxeseni                1        1      1         1    1        0
#> Xyleborus affinis                   0        0      1         0    0        1
#> Xyleborus bispinatus                0        0      0         0    0        0
#> Xyleborus ferrugineus               0        0      0         0    0        1
#> Xyleborus glabratus                 0        0      0         0    0        0
#> Xyleborus posticus                  0        0      0         0    0        0
#> Xyleborus seriatus                  0        0      0         0    1        0
#> Xyleborus spathipennis              0        0      0         0    0        0
#> Xyleborus xylographus               0        0      0         0    0        0
#> Xylosandrus compactus               0        0      0         0    0        0
#> Xylosandrus crassiusculus           0        1      1         0    0        0
#> Xyleborus germanus                  0        0      0         1    1        0
#> Xylosandrus morigerus               0        0      0         0    0        1
#>                           Tetragastris Pouteria Dendropanax Barringtonia
#> Ambrosiodmus obliquus                0        0           0            0
#> Ambrosiodmus rubricollis             0        0           0            0
#> Ambrosiodmus tachygraphus            0        0           0            0
#> Xyleborus atratus                    0        0           0            0
#> Xyleborus obesus                     0        0           0            0
#> mutillatus                           0        0           0            0
#> Coptoborus pseudotenuis              0        0           0            0
#> Euwallacea fornicatus                0        0           0            0
#> Euwallacea validus                   0        0           0            0
#> Sampsonius dampfi                    0        0           0            0
#> Xyleborinus saxeseni                 0        0           0            0
#> Xyleborus affinis                    1        1           1            1
#> Xyleborus bispinatus                 0        0           0            0
#> Xyleborus ferrugineus                0        1           1            1
#> Xyleborus glabratus                  0        0           0            0
#> Xyleborus posticus                   0        0           0            0
#> Xyleborus seriatus                   0        0           0            0
#> Xyleborus spathipennis               0        1           0            0
#> Xyleborus xylographus                0        0           0            0
#> Xylosandrus compactus                0        0           0            0
#> Xylosandrus crassiusculus            0        0           0            0
#> Xyleborus germanus                   0        0           0            0
#> Xylosandrus morigerus                0        1           0            0
#>                           Manilkara Melicoccus Metopium Cajanus Calophyllum
#> Ambrosiodmus obliquus             0          0        0       0           0
#> Ambrosiodmus rubricollis          0          0        0       0           0
#> Ambrosiodmus tachygraphus         0          0        0       0           0
#> Xyleborus atratus                 0          0        0       0           0
#> Xyleborus obesus                  0          0        0       0           0
#> mutillatus                        0          0        0       0           0
#> Coptoborus pseudotenuis           0          0        0       0           0
#> Euwallacea fornicatus             0          0        0       0           0
#> Euwallacea validus                0          0        0       0           0
#> Sampsonius dampfi                 0          0        0       0           0
#> Xyleborinus saxeseni              0          0        0       0           0
#> Xyleborus affinis                 1          1        1       1           1
#> Xyleborus bispinatus              0          0        0       0           0
#> Xyleborus ferrugineus             1          1        0       1           0
#> Xyleborus glabratus               0          0        0       0           0
#> Xyleborus posticus                0          0        0       0           0
#> Xyleborus seriatus                0          0        0       0           0
#> Xyleborus spathipennis            0          0        0       0           0
#> Xyleborus xylographus             0          0        0       0           0
#> Xylosandrus compactus             0          0        0       1           0
#> Xylosandrus crassiusculus         0          0        0       0           1
#> Xyleborus germanus                0          0        0       0           0
#> Xylosandrus morigerus             0          0        0       0           0
#>                           Spathodea Breonia Nectandra Lecythis Hymenaea
#> Ambrosiodmus obliquus             0       0         0        0        0
#> Ambrosiodmus rubricollis          0       0         0        0        0
#> Ambrosiodmus tachygraphus         0       0         0        0        0
#> Xyleborus atratus                 0       0         0        0        0
#> Xyleborus obesus                  0       0         0        0        0
#> mutillatus                        0       0         0        0        0
#> Coptoborus pseudotenuis           0       0         0        0        0
#> Euwallacea fornicatus             0       0         0        0        0
#> Euwallacea validus                0       0         0        0        0
#> Sampsonius dampfi                 0       0         0        0        0
#> Xyleborinus saxeseni              0       0         0        0        0
#> Xyleborus affinis                 1       1         1        1        1
#> Xyleborus bispinatus              0       0         0        0        0
#> Xyleborus ferrugineus             1       0         0        1        1
#> Xyleborus glabratus               0       0         0        0        0
#> Xyleborus posticus                0       0         0        0        0
#> Xyleborus seriatus                0       0         0        0        0
#> Xyleborus spathipennis            0       0         0        0        0
#> Xyleborus xylographus             0       0         0        0        0
#> Xylosandrus compactus             0       0         0        0        0
#> Xylosandrus crassiusculus         0       0         0        0        0
#> Xyleborus germanus                0       0         0        0        0
#> Xylosandrus morigerus             0       0         0        0        0
#>                           Spondias Dacryodes Dracaena Acrocarpus Micropholis
#> Ambrosiodmus obliquus            0         0        0          0           0
#> Ambrosiodmus rubricollis         0         0        0          0           0
#> Ambrosiodmus tachygraphus        0         0        0          0           0
#> Xyleborus atratus                0         0        0          0           0
#> Xyleborus obesus                 0         0        0          0           0
#> mutillatus                       0         0        0          0           0
#> Coptoborus pseudotenuis          0         0        0          0           0
#> Euwallacea fornicatus            0         0        0          0           0
#> Euwallacea validus               0         0        0          0           0
#> Sampsonius dampfi                0         0        0          0           0
#> Xyleborinus saxeseni             0         0        0          0           0
#> Xyleborus affinis                1         1        1          1           1
#> Xyleborus bispinatus             0         0        0          0           0
#> Xyleborus ferrugineus            1         1        1          1           1
#> Xyleborus glabratus              0         0        0          0           0
#> Xyleborus posticus               1         0        1          0           0
#> Xyleborus seriatus               0         0        0          0           0
#> Xyleborus spathipennis           0         0        0          0           0
#> Xyleborus xylographus            0         0        0          0           0
#> Xylosandrus compactus            0         0        0          0           0
#> Xylosandrus crassiusculus        0         0        0          0           0
#> Xyleborus germanus               0         0        0          0           0
#> Xylosandrus morigerus            0         0        0          0           0
#>                           Eschweilera Tabebuia Clethra Alexa Andira Syzygium
#> Ambrosiodmus obliquus               0        0       0     0      0        0
#> Ambrosiodmus rubricollis            0        0       0     0      0        0
#> Ambrosiodmus tachygraphus           0        0       0     0      0        0
#> Xyleborus atratus                   0        0       0     0      0        0
#> Xyleborus obesus                    0        0       0     0      0        0
#> mutillatus                          0        0       0     0      0        0
#> Coptoborus pseudotenuis             0        0       0     0      0        0
#> Euwallacea fornicatus               0        0       0     0      0        0
#> Euwallacea validus                  0        0       0     0      0        0
#> Sampsonius dampfi                   0        0       0     0      0        0
#> Xyleborinus saxeseni                0        0       0     0      0        0
#> Xyleborus affinis                   1        1       1     1      1        1
#> Xyleborus bispinatus                0        0       0     0      0        0
#> Xyleborus ferrugineus               0        1       0     0      1        0
#> Xyleborus glabratus                 0        0       0     0      0        0
#> Xyleborus posticus                  0        0       0     0      0        0
#> Xyleborus seriatus                  0        0       0     0      0        0
#> Xyleborus spathipennis              0        0       0     0      0        0
#> Xyleborus xylographus               0        0       0     0      0        0
#> Xylosandrus compactus               0        0       0     0      0        0
#> Xylosandrus crassiusculus           0        0       0     0      0        1
#> Xyleborus germanus                  0        0       0     0      0        0
#> Xylosandrus morigerus               0        0       0     0      0        0
#>                           Sapium Couma Pentaclethra Cespedesia Sloanea Croton
#> Ambrosiodmus obliquus          0     0            0          0       0      0
#> Ambrosiodmus rubricollis       0     0            0          0       0      0
#> Ambrosiodmus tachygraphus      0     0            0          0       0      0
#> Xyleborus atratus              0     0            0          0       0      0
#> Xyleborus obesus               0     0            0          0       0      0
#> mutillatus                     0     0            0          0       0      0
#> Coptoborus pseudotenuis        0     0            0          0       0      0
#> Euwallacea fornicatus          0     0            0          0       0      0
#> Euwallacea validus             0     0            0          0       0      0
#> Sampsonius dampfi              0     0            0          0       0      0
#> Xyleborinus saxeseni           0     0            0          0       0      0
#> Xyleborus affinis              1     1            1          1       1      1
#> Xyleborus bispinatus           0     0            0          0       0      0
#> Xyleborus ferrugineus          0     1            0          0       0      0
#> Xyleborus glabratus            0     0            0          0       0      0
#> Xyleborus posticus             0     0            0          0       0      0
#> Xyleborus seriatus             0     0            0          0       0      0
#> Xyleborus spathipennis         0     0            0          0       0      0
#> Xyleborus xylographus          0     0            0          0       0      0
#> Xylosandrus compactus          0     0            0          0       0      0
#> Xylosandrus crassiusculus      0     0            0          0       0      0
#> Xyleborus germanus             0     0            0          0       0      0
#> Xylosandrus morigerus          0     0            0          0       0      0
#>                           Cocos Saccharum Cecropia Humiriastrum Trichilia
#> Ambrosiodmus obliquus         0         0        0            0         0
#> Ambrosiodmus rubricollis      0         0        0            0         0
#> Ambrosiodmus tachygraphus     0         0        0            0         0
#> Xyleborus atratus             0         0        0            0         0
#> Xyleborus obesus              0         0        0            0         0
#> mutillatus                    0         0        0            0         0
#> Coptoborus pseudotenuis       0         0        0            0         0
#> Euwallacea fornicatus         0         0        0            0         0
#> Euwallacea validus            0         0        0            0         0
#> Sampsonius dampfi             0         0        0            0         0
#> Xyleborinus saxeseni          0         0        0            0         0
#> Xyleborus affinis             1         1        1            1         1
#> Xyleborus bispinatus          0         0        0            0         0
#> Xyleborus ferrugineus         0         0        1            1         0
#> Xyleborus glabratus           0         0        0            0         0
#> Xyleborus posticus            0         0        1            0         0
#> Xyleborus seriatus            0         0        0            0         0
#> Xyleborus spathipennis        0         0        0            0         0
#> Xyleborus xylographus         0         0        0            0         0
#> Xylosandrus compactus         0         0        0            0         0
#> Xylosandrus crassiusculus     0         1        1            0         0
#> Xyleborus germanus            0         0        0            0         0
#> Xylosandrus morigerus         0         0        1            0         0
#>                           Leucaena Toulicia Pachira Cyrilla Syagrus Liquidambar
#> Ambrosiodmus obliquus            0        0       0       0       0           0
#> Ambrosiodmus rubricollis         0        0       0       0       0           0
#> Ambrosiodmus tachygraphus        0        0       0       0       0           0
#> Xyleborus atratus                0        0       0       0       0           0
#> Xyleborus obesus                 0        0       0       0       0           0
#> mutillatus                       0        0       0       0       0           0
#> Coptoborus pseudotenuis          0        0       0       0       0           0
#> Euwallacea fornicatus            0        0       0       0       0           0
#> Euwallacea validus               0        0       0       0       0           0
#> Sampsonius dampfi                0        0       0       0       0           0
#> Xyleborinus saxeseni             0        0       0       0       0           0
#> Xyleborus affinis                1        1       1       1       1           1
#> Xyleborus bispinatus             0        0       0       0       0           0
#> Xyleborus ferrugineus            0        0       1       0       0           0
#> Xyleborus glabratus              0        0       0       0       0           0
#> Xyleborus posticus               1        0       0       1       0           0
#> Xyleborus seriatus               0        0       0       0       0           0
#> Xyleborus spathipennis           0        0       0       0       0           0
#> Xyleborus xylographus            0        0       0       0       0           0
#> Xylosandrus compactus            0        0       0       0       0           1
#> Xylosandrus crassiusculus        0        0       0       0       0           1
#> Xyleborus germanus               0        0       0       0       0           0
#> Xylosandrus morigerus            0        0       0       0       0           0
#>                           Buchenavia Inga Thouinidium Aralia Fissicalyx
#> Ambrosiodmus obliquus              0    0           0      0          0
#> Ambrosiodmus rubricollis           0    0           0      0          0
#> Ambrosiodmus tachygraphus          0    0           0      0          0
#> Xyleborus atratus                  0    0           0      0          0
#> Xyleborus obesus                   0    0           0      0          0
#> mutillatus                         0    0           0      0          0
#> Coptoborus pseudotenuis            0    0           0      0          0
#> Euwallacea fornicatus              0    0           0      0          0
#> Euwallacea validus                 0    0           0      0          0
#> Sampsonius dampfi                  0    0           0      0          0
#> Xyleborinus saxeseni               0    0           0      0          0
#> Xyleborus affinis                  1    1           0      0          0
#> Xyleborus bispinatus               0    0           0      0          0
#> Xyleborus ferrugineus              0    0           1      1          1
#> Xyleborus glabratus                0    0           0      0          0
#> Xyleborus posticus                 0    0           0      0          0
#> Xyleborus seriatus                 0    0           0      0          0
#> Xyleborus spathipennis             0    0           0      0          0
#> Xyleborus xylographus              0    0           0      0          0
#> Xylosandrus compactus              0    0           0      0          0
#> Xylosandrus crassiusculus          0    0           0      0          0
#> Xyleborus germanus                 0    0           0      0          0
#> Xylosandrus morigerus              0    0           0      0          0
#>                           Avicennia Vitex Aspidosperma Lonchocarpus Picea
#> Ambrosiodmus obliquus             0     0            0            0     0
#> Ambrosiodmus rubricollis          0     0            0            0     0
#> Ambrosiodmus tachygraphus         0     0            0            0     0
#> Xyleborus atratus                 0     0            0            0     0
#> Xyleborus obesus                  0     0            0            0     0
#> mutillatus                        0     0            0            0     0
#> Coptoborus pseudotenuis           0     0            0            0     0
#> Euwallacea fornicatus             0     0            0            0     0
#> Euwallacea validus                0     0            0            0     0
#> Sampsonius dampfi                 0     0            0            0     0
#> Xyleborinus saxeseni              0     0            0            0     0
#> Xyleborus affinis                 0     0            0            0     0
#> Xyleborus bispinatus              0     0            0            0     0
#> Xyleborus ferrugineus             1     1            1            1     0
#> Xyleborus glabratus               0     0            0            0     1
#> Xyleborus posticus                0     0            0            0     0
#> Xyleborus seriatus                0     0            0            0     0
#> Xyleborus spathipennis            0     0            0            0     0
#> Xyleborus xylographus             0     0            0            0     0
#> Xylosandrus compactus             0     0            0            0     0
#> Xylosandrus crassiusculus         0     1            0            0     0
#> Xyleborus germanus                0     0            0            0     0
#> Xylosandrus morigerus             0     0            0            0     0
#>                           Litsea Phoebe Schima Astronium Larix Kalopanax Thuja
#> Ambrosiodmus obliquus          0      0      0         0     0         0     0
#> Ambrosiodmus rubricollis       0      0      0         0     0         0     0
#> Ambrosiodmus tachygraphus      0      0      0         0     0         0     0
#> Xyleborus atratus              0      0      0         0     0         0     0
#> Xyleborus obesus               0      0      0         0     0         0     0
#> mutillatus                     0      0      0         0     0         0     0
#> Coptoborus pseudotenuis        0      0      0         0     0         0     0
#> Euwallacea fornicatus          0      0      0         0     0         0     0
#> Euwallacea validus             0      0      0         0     0         0     0
#> Sampsonius dampfi              0      0      0         0     0         0     0
#> Xyleborinus saxeseni           0      0      0         0     0         0     0
#> Xyleborus affinis              0      0      0         0     0         0     0
#> Xyleborus bispinatus           0      0      0         0     0         0     0
#> Xyleborus ferrugineus          0      0      0         0     0         0     0
#> Xyleborus glabratus            1      1      1         0     0         0     0
#> Xyleborus posticus             0      0      0         1     0         0     0
#> Xyleborus seriatus             0      0      0         0     1         1     1
#> Xyleborus spathipennis         0      0      0         0     0         0     0
#> Xyleborus xylographus          0      0      0         0     0         0     0
#> Xylosandrus compactus          0      0      0         0     0         0     0
#> Xylosandrus crassiusculus      0      1      0         0     0         0     0
#> Xyleborus germanus             0      0      1         0     0         0     0
#> Xylosandrus morigerus          0      0      0         1     0         0     0
#>                           Aesculus Ilex Gleditsia Coffea Arundina Myrica
#> Ambrosiodmus obliquus            0    0         0      0        0      0
#> Ambrosiodmus rubricollis         0    0         0      0        0      0
#> Ambrosiodmus tachygraphus        0    0         0      0        0      0
#> Xyleborus atratus                0    0         0      0        0      0
#> Xyleborus obesus                 0    0         0      0        0      0
#> mutillatus                       0    0         0      0        0      0
#> Coptoborus pseudotenuis          0    0         0      0        0      0
#> Euwallacea fornicatus            0    0         0      0        0      0
#> Euwallacea validus               0    0         0      0        0      0
#> Sampsonius dampfi                0    0         0      0        0      0
#> Xyleborinus saxeseni             0    0         0      0        0      0
#> Xyleborus affinis                0    0         0      0        0      0
#> Xyleborus bispinatus             0    0         0      0        0      0
#> Xyleborus ferrugineus            0    0         0      0        0      0
#> Xyleborus glabratus              0    0         0      0        0      0
#> Xyleborus posticus               0    0         0      0        0      0
#> Xyleborus seriatus               1    0         0      0        0      0
#> Xyleborus spathipennis           0    0         0      0        0      0
#> Xyleborus xylographus            0    0         0      0        0      0
#> Xylosandrus compactus            0    1         1      1        1      1
#> Xylosandrus crassiusculus        0    0         1      0        0      0
#> Xyleborus germanus               1    0         0      0        0      0
#> Xylosandrus morigerus            0    0         0      1        0      0
#>                           Piscidia Olea Piper Ardisia Koelreuteria Zanthoxylum
#> Ambrosiodmus obliquus            0    0     0       0            0           0
#> Ambrosiodmus rubricollis         0    0     0       0            0           0
#> Ambrosiodmus tachygraphus        0    0     0       0            0           0
#> Xyleborus atratus                0    0     0       0            0           0
#> Xyleborus obesus                 0    0     0       0            0           0
#> mutillatus                       0    0     0       0            0           0
#> Coptoborus pseudotenuis          0    0     0       0            0           0
#> Euwallacea fornicatus            0    0     0       0            0           0
#> Euwallacea validus               0    0     0       0            0           0
#> Sampsonius dampfi                0    0     0       0            0           0
#> Xyleborinus saxeseni             0    0     0       0            0           0
#> Xyleborus affinis                0    0     0       0            0           0
#> Xyleborus bispinatus             0    0     0       0            0           0
#> Xyleborus ferrugineus            0    0     0       0            0           0
#> Xyleborus glabratus              0    0     0       0            0           0
#> Xyleborus posticus               0    0     0       0            0           0
#> Xyleborus seriatus               0    0     0       0            0           0
#> Xyleborus spathipennis           0    0     0       0            0           0
#> Xyleborus xylographus            0    0     0       0            0           0
#> Xylosandrus compactus            1    1     1       1            1           0
#> Xylosandrus crassiusculus        0    0     0       0            1           1
#> Xyleborus germanus               0    0     0       0            0           0
#> Xylosandrus morigerus            0    0     0       0            0           0
#>                           Pycnanthus Holigarna Dypsis Ipomoea Hopea Styrax
#> Ambrosiodmus obliquus              0         0      0       0     0      0
#> Ambrosiodmus rubricollis           0         0      0       0     0      0
#> Ambrosiodmus tachygraphus          0         0      0       0     0      0
#> Xyleborus atratus                  0         0      0       0     0      0
#> Xyleborus obesus                   0         0      0       0     0      0
#> mutillatus                         0         0      0       0     0      0
#> Coptoborus pseudotenuis            0         0      0       0     0      0
#> Euwallacea fornicatus              0         0      0       0     0      0
#> Euwallacea validus                 0         0      0       0     0      0
#> Sampsonius dampfi                  0         0      0       0     0      0
#> Xyleborinus saxeseni               0         0      0       0     0      0
#> Xyleborus affinis                  0         0      0       0     0      0
#> Xyleborus bispinatus               0         0      0       0     0      0
#> Xyleborus ferrugineus              0         0      0       0     0      0
#> Xyleborus glabratus                0         0      0       0     0      0
#> Xyleborus posticus                 0         0      0       0     0      0
#> Xyleborus seriatus                 0         0      0       0     0      0
#> Xyleborus spathipennis             0         0      0       0     0      0
#> Xyleborus xylographus              0         0      0       0     0      0
#> Xylosandrus compactus              0         0      0       0     0      0
#> Xylosandrus crassiusculus          1         1      1       1     1      1
#> Xyleborus germanus                 0         0      0       0     0      1
#> Xylosandrus morigerus              0         0      0       0     0      0
#>                           Pourouma Sambucus Firmiana Streblus Castilla
#> Ambrosiodmus obliquus            0        0        0        0        0
#> Ambrosiodmus rubricollis         0        0        0        0        0
#> Ambrosiodmus tachygraphus        0        0        0        0        0
#> Xyleborus atratus                0        0        0        0        0
#> Xyleborus obesus                 0        0        0        0        0
#> mutillatus                       0        0        0        0        0
#> Coptoborus pseudotenuis          0        0        0        0        0
#> Euwallacea fornicatus            0        0        0        0        0
#> Euwallacea validus               0        0        0        0        0
#> Sampsonius dampfi                0        0        0        0        0
#> Xyleborinus saxeseni             0        0        0        0        0
#> Xyleborus affinis                0        0        0        0        0
#> Xyleborus bispinatus             0        0        0        0        0
#> Xyleborus ferrugineus            0        0        0        0        0
#> Xyleborus glabratus              0        0        0        0        0
#> Xyleborus posticus               0        0        0        0        0
#> Xyleborus seriatus               0        0        0        0        0
#> Xyleborus spathipennis           0        0        0        0        0
#> Xyleborus xylographus            0        0        0        0        0
#> Xylosandrus compactus            0        0        0        0        0
#> Xylosandrus crassiusculus        1        1        1        1        1
#> Xyleborus germanus               0        0        0        0        0
#> Xylosandrus morigerus            0        1        0        0        0
#>                           Castanopsis Vochysia Mesua Ongokea Tectona
#> Ambrosiodmus obliquus               0        0     0       0       0
#> Ambrosiodmus rubricollis            0        0     0       0       0
#> Ambrosiodmus tachygraphus           0        0     0       0       0
#> Xyleborus atratus                   0        0     0       0       0
#> Xyleborus obesus                    0        0     0       0       0
#> mutillatus                          0        0     0       0       0
#> Coptoborus pseudotenuis             0        0     0       0       0
#> Euwallacea fornicatus               0        0     0       0       0
#> Euwallacea validus                  0        0     0       0       0
#> Sampsonius dampfi                   0        0     0       0       0
#> Xyleborinus saxeseni                0        0     0       0       0
#> Xyleborus affinis                   0        0     0       0       0
#> Xyleborus bispinatus                0        0     0       0       0
#> Xyleborus ferrugineus               0        0     0       0       0
#> Xyleborus glabratus                 0        0     0       0       0
#> Xyleborus posticus                  0        0     0       0       0
#> Xyleborus seriatus                  0        0     0       0       0
#> Xyleborus spathipennis              0        0     0       0       0
#> Xyleborus xylographus               0        0     0       0       0
#> Xylosandrus compactus               0        0     0       0       0
#> Xylosandrus crassiusculus           1        1     1       1       1
#> Xyleborus germanus                  0        0     0       0       0
#> Xylosandrus morigerus               1        0     0       0       1
#>                           Erythrophleum Vateria Khaya Bischofia Aucoumea
#> Ambrosiodmus obliquus                 0       0     0         0        0
#> Ambrosiodmus rubricollis              0       0     0         0        0
#> Ambrosiodmus tachygraphus             0       0     0         0        0
#> Xyleborus atratus                     0       0     0         0        0
#> Xyleborus obesus                      0       0     0         0        0
#> mutillatus                            0       0     0         0        0
#> Coptoborus pseudotenuis               0       0     0         0        0
#> Euwallacea fornicatus                 0       0     0         0        0
#> Euwallacea validus                    0       0     0         0        0
#> Sampsonius dampfi                     0       0     0         0        0
#> Xyleborinus saxeseni                  0       0     0         0        0
#> Xyleborus affinis                     0       0     0         0        0
#> Xyleborus bispinatus                  0       0     0         0        0
#> Xyleborus ferrugineus                 0       0     0         0        0
#> Xyleborus glabratus                   0       0     0         0        0
#> Xyleborus posticus                    0       0     0         0        0
#> Xyleborus seriatus                    0       0     0         0        0
#> Xyleborus spathipennis                0       0     0         0        0
#> Xyleborus xylographus                 0       0     0         0        0
#> Xylosandrus compactus                 0       0     0         0        0
#> Xylosandrus crassiusculus             1       1     1         1        1
#> Xyleborus germanus                    0       0     0         0        0
#> Xylosandrus morigerus                 0       0     0         0        0
#>                           Murraya Leplaea Sageraea Blakea Dillenia Leea
#> Ambrosiodmus obliquus           0       0        0      0        0    0
#> Ambrosiodmus rubricollis        0       0        0      0        0    0
#> Ambrosiodmus tachygraphus       0       0        0      0        0    0
#> Xyleborus atratus               0       0        0      0        0    0
#> Xyleborus obesus                0       0        0      0        0    0
#> mutillatus                      0       0        0      0        0    0
#> Coptoborus pseudotenuis         0       0        0      0        0    0
#> Euwallacea fornicatus           0       0        0      0        0    0
#> Euwallacea validus              0       0        0      0        0    0
#> Sampsonius dampfi               0       0        0      0        0    0
#> Xyleborinus saxeseni            0       0        0      0        0    0
#> Xyleborus affinis               0       0        0      0        0    0
#> Xyleborus bispinatus            0       0        0      0        0    0
#> Xyleborus ferrugineus           0       0        0      0        0    0
#> Xyleborus glabratus             0       0        0      0        0    0
#> Xyleborus posticus              0       0        0      0        0    0
#> Xyleborus seriatus              0       0        0      0        0    0
#> Xyleborus spathipennis          0       0        0      0        0    0
#> Xyleborus xylographus           0       0        0      0        0    0
#> Xylosandrus compactus           0       0        0      0        0    0
#> Xylosandrus crassiusculus       1       1        1      1        1    1
#> Xyleborus germanus              0       0        0      0        0    0
#> Xylosandrus morigerus           0       0        0      0        0    0
#>                           Cannabis Elaeocarpus Swintonia Lagerstroemia Ungnadia
#> Ambrosiodmus obliquus            0           0         0             0        0
#> Ambrosiodmus rubricollis         0           0         0             0        0
#> Ambrosiodmus tachygraphus        0           0         0             0        0
#> Xyleborus atratus                0           0         0             0        0
#> Xyleborus obesus                 0           0         0             0        0
#> mutillatus                       0           0         0             0        0
#> Coptoborus pseudotenuis          0           0         0             0        0
#> Euwallacea fornicatus            0           0         0             0        0
#> Euwallacea validus               0           0         0             0        0
#> Sampsonius dampfi                0           0         0             0        0
#> Xyleborinus saxeseni             0           0         0             0        0
#> Xyleborus affinis                0           0         0             0        0
#> Xyleborus bispinatus             0           0         0             0        0
#> Xyleborus ferrugineus            0           0         0             0        0
#> Xyleborus glabratus              0           0         0             0        0
#> Xyleborus posticus               0           0         0             0        0
#> Xyleborus seriatus               0           0         0             0        0
#> Xyleborus spathipennis           0           0         0             0        0
#> Xyleborus xylographus            0           0         0             0        0
#> Xylosandrus compactus            0           0         0             0        0
#> Xylosandrus crassiusculus        1           1         1             1        1
#> Xyleborus germanus               0           0         0             0        0
#> Xylosandrus morigerus            0           0         0             0        0
#>                           Chloroxylon Wrightia Toona Gluta Sterculia
#> Ambrosiodmus obliquus               0        0     0     0         0
#> Ambrosiodmus rubricollis            0        0     0     0         0
#> Ambrosiodmus tachygraphus           0        0     0     0         0
#> Xyleborus atratus                   0        0     0     0         0
#> Xyleborus obesus                    0        0     0     0         0
#> mutillatus                          0        0     0     0         0
#> Coptoborus pseudotenuis             0        0     0     0         0
#> Euwallacea fornicatus               0        0     0     0         0
#> Euwallacea validus                  0        0     0     0         0
#> Sampsonius dampfi                   0        0     0     0         0
#> Xyleborinus saxeseni                0        0     0     0         0
#> Xyleborus affinis                   0        0     0     0         0
#> Xyleborus bispinatus                0        0     0     0         0
#> Xyleborus ferrugineus               0        0     0     0         0
#> Xyleborus glabratus                 0        0     0     0         0
#> Xyleborus posticus                  0        0     0     0         0
#> Xyleborus seriatus                  0        0     0     0         0
#> Xyleborus spathipennis              0        0     0     0         0
#> Xyleborus xylographus               0        0     0     0         0
#> Xylosandrus compactus               0        0     0     0         0
#> Xylosandrus crassiusculus           1        1     1     1         1
#> Xyleborus germanus                  0        0     0     0         0
#> Xylosandrus morigerus               0        0     0     0         0
#>                           Lithocarpus Scorodophloeus Dipterocarpus Sassafras
#> Ambrosiodmus obliquus               0              0             0         0
#> Ambrosiodmus rubricollis            0              0             0         0
#> Ambrosiodmus tachygraphus           0              0             0         0
#> Xyleborus atratus                   0              0             0         0
#> Xyleborus obesus                    0              0             0         0
#> mutillatus                          0              0             0         0
#> Coptoborus pseudotenuis             0              0             0         0
#> Euwallacea fornicatus               0              0             0         0
#> Euwallacea validus                  0              0             0         0
#> Sampsonius dampfi                   0              0             0         0
#> Xyleborinus saxeseni                0              0             0         0
#> Xyleborus affinis                   0              0             0         0
#> Xyleborus bispinatus                0              0             0         0
#> Xyleborus ferrugineus               0              0             0         0
#> Xyleborus glabratus                 0              0             0         0
#> Xyleborus posticus                  0              0             0         0
#> Xyleborus seriatus                  0              0             0         0
#> Xyleborus spathipennis              0              0             0         0
#> Xyleborus xylographus               0              0             0         0
#> Xylosandrus compactus               0              0             0         0
#> Xylosandrus crassiusculus           1              1             1         0
#> Xyleborus germanus                  0              0             0         1
#> Xylosandrus morigerus               0              0             0         0
#>                           Fraxinus Nyssa Salix Ziziphus Platanus Juniperus Cola
#> Ambrosiodmus obliquus            0     0     0        0        0         0    0
#> Ambrosiodmus rubricollis         0     0     0        0        0         0    0
#> Ambrosiodmus tachygraphus        0     0     0        0        0         0    0
#> Xyleborus atratus                0     0     0        0        0         0    0
#> Xyleborus obesus                 0     0     0        0        0         0    0
#> mutillatus                       0     0     0        0        0         0    0
#> Coptoborus pseudotenuis          0     0     0        0        0         0    0
#> Euwallacea fornicatus            0     0     0        0        0         0    0
#> Euwallacea validus               0     0     0        0        0         0    0
#> Sampsonius dampfi                0     0     0        0        0         0    0
#> Xyleborinus saxeseni             0     0     0        0        0         0    0
#> Xyleborus affinis                0     0     0        0        0         0    0
#> Xyleborus bispinatus             0     0     0        0        0         0    0
#> Xyleborus ferrugineus            0     0     0        0        0         0    0
#> Xyleborus glabratus              0     0     0        0        0         0    0
#> Xyleborus posticus               0     0     0        0        0         0    0
#> Xyleborus seriatus               0     0     0        0        0         0    0
#> Xyleborus spathipennis           0     0     0        0        0         0    0
#> Xyleborus xylographus            0     0     0        0        0         0    0
#> Xylosandrus compactus            0     0     0        0        0         0    0
#> Xylosandrus crassiusculus        0     0     0        0        0         0    0
#> Xyleborus germanus               1     1     1        1        1         1    0
#> Xylosandrus morigerus            0     0     0        0        0         0    1
#>                           Machaerium Altingia Gymnopodium Schoutenia Acaciella
#> Ambrosiodmus obliquus              0        0           0          0         0
#> Ambrosiodmus rubricollis           0        0           0          0         0
#> Ambrosiodmus tachygraphus          0        0           0          0         0
#> Xyleborus atratus                  0        0           0          0         0
#> Xyleborus obesus                   0        0           0          0         0
#> mutillatus                         0        0           0          0         0
#> Coptoborus pseudotenuis            0        0           0          0         0
#> Euwallacea fornicatus              0        0           0          0         0
#> Euwallacea validus                 0        0           0          0         0
#> Sampsonius dampfi                  0        0           0          0         0
#> Xyleborinus saxeseni               0        0           0          0         0
#> Xyleborus affinis                  0        0           0          0         0
#> Xyleborus bispinatus               0        0           0          0         0
#> Xyleborus ferrugineus              0        0           0          0         0
#> Xyleborus glabratus                0        0           0          0         0
#> Xyleborus posticus                 0        0           0          0         0
#> Xyleborus seriatus                 0        0           0          0         0
#> Xyleborus spathipennis             0        0           0          0         0
#> Xyleborus xylographus              0        0           0          0         0
#> Xylosandrus compactus              0        0           0          0         0
#> Xylosandrus crassiusculus          0        0           0          0         0
#> Xyleborus germanus                 0        0           0          0         0
#> Xylosandrus morigerus              1        1           1          1         1
#>                           Licania Tarennoidea Grewia Derris Butea Macrolenes
#> Ambrosiodmus obliquus           0           0      0      0     0          0
#> Ambrosiodmus rubricollis        0           0      0      0     0          0
#> Ambrosiodmus tachygraphus       0           0      0      0     0          0
#> Xyleborus atratus               0           0      0      0     0          0
#> Xyleborus obesus                0           0      0      0     0          0
#> mutillatus                      0           0      0      0     0          0
#> Coptoborus pseudotenuis         0           0      0      0     0          0
#> Euwallacea fornicatus           0           0      0      0     0          0
#> Euwallacea validus              0           0      0      0     0          0
#> Sampsonius dampfi               0           0      0      0     0          0
#> Xyleborinus saxeseni            0           0      0      0     0          0
#> Xyleborus affinis               0           0      0      0     0          0
#> Xyleborus bispinatus            0           0      0      0     0          0
#> Xyleborus ferrugineus           0           0      0      0     0          0
#> Xyleborus glabratus             0           0      0      0     0          0
#> Xyleborus posticus              0           0      0      0     0          0
#> Xyleborus seriatus              0           0      0      0     0          0
#> Xyleborus spathipennis          0           0      0      0     0          0
#> Xyleborus xylographus           0           0      0      0     0          0
#> Xylosandrus compactus           0           0      0      0     0          0
#> Xylosandrus crassiusculus       0           0      0      0     0          0
#> Xyleborus germanus              0           0      0      0     0          0
#> Xylosandrus morigerus           1           1      1      1     1          1
#>                           Schleichera Bixa Austroeupatorium Schizolobium
#> Ambrosiodmus obliquus               0    0                0            0
#> Ambrosiodmus rubricollis            0    0                0            0
#> Ambrosiodmus tachygraphus           0    0                0            0
#> Xyleborus atratus                   0    0                0            0
#> Xyleborus obesus                    0    0                0            0
#> mutillatus                          0    0                0            0
#> Coptoborus pseudotenuis             0    0                0            0
#> Euwallacea fornicatus               0    0                0            0
#> Euwallacea validus                  0    0                0            0
#> Sampsonius dampfi                   0    0                0            0
#> Xyleborinus saxeseni                0    0                0            0
#> Xyleborus affinis                   0    0                0            0
#> Xyleborus bispinatus                0    0                0            0
#> Xyleborus ferrugineus               0    0                0            0
#> Xyleborus glabratus                 0    0                0            0
#> Xyleborus posticus                  0    0                0            0
#> Xyleborus seriatus                  0    0                0            0
#> Xyleborus spathipennis              0    0                0            0
#> Xyleborus xylographus               0    0                0            0
#> Xylosandrus compactus               0    0                0            0
#> Xylosandrus crassiusculus           0    0                0            0
#> Xyleborus germanus                  0    0                0            0
#> Xylosandrus morigerus               1    1                1            1
#>                           Adenanthera Ceiba Esenbeckia Dendrobium Centrosema
#> Ambrosiodmus obliquus               0     0          0          0          0
#> Ambrosiodmus rubricollis            0     0          0          0          0
#> Ambrosiodmus tachygraphus           0     0          0          0          0
#> Xyleborus atratus                   0     0          0          0          0
#> Xyleborus obesus                    0     0          0          0          0
#> mutillatus                          0     0          0          0          0
#> Coptoborus pseudotenuis             0     0          0          0          0
#> Euwallacea fornicatus               0     0          0          0          0
#> Euwallacea validus                  0     0          0          0          0
#> Sampsonius dampfi                   0     0          0          0          0
#> Xyleborinus saxeseni                0     0          0          0          0
#> Xyleborus affinis                   0     0          0          0          0
#> Xyleborus bispinatus                0     0          0          0          0
#> Xyleborus ferrugineus               0     0          0          0          0
#> Xyleborus glabratus                 0     0          0          0          0
#> Xyleborus posticus                  0     0          0          0          0
#> Xyleborus seriatus                  0     0          0          0          0
#> Xyleborus spathipennis              0     0          0          0          0
#> Xyleborus xylographus               0     0          0          0          0
#> Xylosandrus compactus               0     0          0          0          0
#> Xylosandrus crassiusculus           0     0          0          0          0
#> Xyleborus germanus                  0     0          0          0          0
#> Xylosandrus morigerus               1     1          1          1          1
#>                           Claoxylon Eugenia Flemingia Alseis Eusideroxylon
#> Ambrosiodmus obliquus             0       0         0      0             0
#> Ambrosiodmus rubricollis          0       0         0      0             0
#> Ambrosiodmus tachygraphus         0       0         0      0             0
#> Xyleborus atratus                 0       0         0      0             0
#> Xyleborus obesus                  0       0         0      0             0
#> mutillatus                        0       0         0      0             0
#> Coptoborus pseudotenuis           0       0         0      0             0
#> Euwallacea fornicatus             0       0         0      0             0
#> Euwallacea validus                0       0         0      0             0
#> Sampsonius dampfi                 0       0         0      0             0
#> Xyleborinus saxeseni              0       0         0      0             0
#> Xyleborus affinis                 0       0         0      0             0
#> Xyleborus bispinatus              0       0         0      0             0
#> Xyleborus ferrugineus             0       0         0      0             0
#> Xyleborus glabratus               0       0         0      0             0
#> Xyleborus posticus                0       0         0      0             0
#> Xyleborus seriatus                0       0         0      0             0
#> Xyleborus spathipennis            0       0         0      0             0
#> Xyleborus xylographus             0       0         0      0             0
#> Xylosandrus compactus             0       0         0      0             0
#> Xylosandrus crassiusculus         0       0         0      0             0
#> Xyleborus germanus                0       0         0      0             0
#> Xylosandrus morigerus             1       1         1      1             1
```
