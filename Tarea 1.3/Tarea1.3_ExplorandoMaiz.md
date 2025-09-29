# Preguntas de la Tarea 1.3

## Pregunta 1: ¿Qué tipo de objeto creamos al cargar la base?

    maiz <- read.delim("C:/Users/pame5/Documents/BioinfinvRepro/Unidad1/Sesion3/PracUni1Ses3/maices/meta/maizteocintle_SNP50k_meta_extended.txt", fileEncoding = "latin1")
    class(maiz)
    
    ## [1] "data.frame"

## Pregunta 2 ¿Cómo se ven las primeras 6 líneas del archivo?

    head(maiz) 
    
    ##   OrderColecta NSiembra          Origen                  Raza     Estado
    ## 1            2        3  INIFAP-2009-16              Apachito  Chihuahua
    ## 2            1       68         2009-72            ChalqueÌ±o   Tlaxcala
    ## 3            1       91         2007-33 Dulcillo del Noroeste     Sonora
    ## 4            2       39    Chis-2009-18            Dzit-Bacal    Chiapas
    ## 5            2       12 Celaya-2009-114                Celaya Guanajuato
    ## 6            2       41   Celaya-2009-2   Elotes Occidentales Guanajuato
    ##     Num_Colecta       Nombre_comun         Raza_Primaria Raza_Secundaria
    ## 1            16      Ocho carreras              Apachito              NA
    ## 2            22 ChalqueÌ±o Criollo            ChalqueÌ±o              NA
    ## 3   SON2007-033        MaÌ_z Dulce Dulcillo del Noroeste              NA
    ## 4   Repetida 18     Olotillo crema            Dzit-Bacal              NA
    ## 5 2009-REPO-114      MaÌ_z Criollo                Celaya              NA
    ## 6 2009-REPO-002      MaÌ_z Criollo   Elotes Occidentales              NA
    ##   AÌ.o._de_colecta              Localidad  Municipio   Estado.1   Longitud
    ## 1             2009           Santo TomÌÁs   Guerrero  Chihuahua -107.58300
    ## 2             2008       Ignacio Zaragoza Cuapiaxtla   Tlaxcala  -97.92944
    ## 3             2007            Agua Blanca    YÌ©cora     Sonora -108.92467
    ## 4             2009 Nuevo Vicente Guerrero Villacorzo    Chiapas  -92.97972
    ## 5             2009            El Ahuacate  Uriangato Guanajuato -101.11889
    ## 6             2009              Comonfort  Comonfort Guanajuato -100.76583
    ##    Latitud Altitud                         Ruizetal2008_grupo
    ## 1 28.68578    1975                       1A_templado540-640mm
    ## 2 19.29083    2497                                           
    ## 3 28.53703    1435                     2A_semicalido500-870mm
    ## 4 16.03225     618                     3A_muycalido990-1360mm
    ## 5 20.09111    1877 2A_semicalido500-870mm y 1B_templado>650mm
    ## 6 20.74417    1800 2C_semicalido740-855mm y 1B_templado>650mm
    ##     Sanchezetal_grupo Categ.Altitud                 Rbiogeo
    ## 1 Sierra de Chihuahua           mid Sierra Madre Occidental
    ## 2             CÌ_nico          high           Eje Volcanico
    ## 3           Chapalote           mid Sierra Madre Occidental
    ## 4 MaduraciÌ_n tardÌ_a           low      Costa del Pacifico
    ## 5 Dentados tropicales           mid           Eje Volcanico
    ## 6        Ocho Hileras           mid           Eje Volcanico
    ##              DivFloristic  PeralesBiog
    ## 1 Sierra Madre Occidental Caå_ones Chi
    ## 2  Serranias Meridionales  Mesa Centra
    ## 3 Sierra Madre Occidental  Sierras del
    ## 4  Serranias Transismicas      Chiapas
    ## 5            Altiplanicie       BajåÁo
    ## 6            Altiplanicie       BajåÁo

## Pregunta 3 ¿Cuántas muestras hay?

    nrow(maiz) 
    
    ## [1] 165

## Pregunta 4 ¿De cuántos estados se tienen muestras?

    length(unique(maiz$Estado)) 
    
    ## [1] 19

## Pregunta 5 ¿Cuántas muestras fueron colectadas antes de 1980?

    sum(maiz$A_o._de_colecta < 1980, na.rm = TRUE) 
    
    ## [1] 0

## Pregunta 6 ¿Cuántas muestras hay de cada raza?

    table(maiz$Raza)
    
    ## 
    ##                       Ancho                    Apachito 
    ##                           3                           2 
    ##                   Arrocillo                        Azul 
    ##                           4                           2 
    ##            Blando de Sonora                        Bofo 
    ##                           1                           1 
    ##               Cacahuacintle                      Celaya 
    ##                           5                           3 
    ##                  ChalqueÌ±o                   Chapalote 
    ##                           7                           2 
    ##                     CÌ_nico            CÌ_nico NorteÌ±o 
    ##                          16                           3 
    ##                    Comiteco Complejo Serrano de Jalisco 
    ##                           5                           2 
    ##                      Conejo                Coscomatepec 
    ##                           4                           3 
    ##     Cristalino de Chihuahua                       Dulce 
    ##                           2                           1 
    ##       Dulcillo del Noroeste                  Dzit-Bacal 
    ##                           2                           3 
    ##          Elotero de Sinaloa             Elotes CÌ_nicos 
    ##                           5                          14 
    ##         Elotes Occidentales                       Gordo 
    ##                           4                           2 
    ##                        Jala                     Mushito 
    ##                           4                           3 
    ##           Nal-tel de Altura                     OlotÌ_n 
    ##                           5                           4 
    ##                    Olotillo                    OnaveÌ±o 
    ##                           6                           2 
    ##       Palomero de Chihuahua         Palomero ToluqueÌ±o 
    ##                           1                           1 
    ##                   Pepitilla                      RatÌ_n 
    ##                           4                           3 
    ##                  Reventador            Tablilla de Ocho 
    ##                           2                           2 
    ##                 Tabloncillo           Tabloncillo Perla 
    ##                           4                           3 
    ##                       Tehua                  Tepecintle 
    ##                           2                           4 
    ##                    TuxpeÌ±o           TuxpeÌ±o NorteÌ±o 
    ##                           4                           2 
    ##                    VandeÌ±o           Zamorano Amarillo 
    ##                           4                           3 
    ##              Zapalote Chico             Zapalote Grande 
    ##                           1                           1 
    ##             Zea m. mexicana          Zea m. parviglumis 
    ##                           2                           2

## Pregunta 7 En promedio, ¿a qué altitud fueron colectadas las muestras?

    mean(maiz$Altitud, na.rm = TRUE)
    
    ## [1] 1519.242

## Pregunta 8 ¿A qué altitud máxima y mínima fueron colectadas?

    max(maiz$Altitud, na.rm = TRUE)
    
    ## [1] 2769
    
    min(maiz$Altitud, na.rm = TRUE)
    
    ## [1] 5

## Pregunta 9 Crear un nuevo data frame solo con las muestras de la raza Olotillo

    olotillo <- subset(maiz, Raza == "Olotillo")
    head(olotillo)
    
    ##     OrderColecta NSiembra   Origen     Raza  Estado Num_Colecta
    ## 68             1       76 2009-136 Olotillo Hidalgo         136
    ## 72             1       77 2009-152 Olotillo Hidalgo         152
    ## 75             1       78  2009-18 Olotillo Chiapas          18
    ## 79             1       79  2009-31 Olotillo Chiapas          31
    ## 123            1      165 2009-260 Olotillo Chiapas         260
    ## 159            1      200 2009-264 Olotillo Chiapas         264
    ##                      Nombre_comun Raza_Primaria Raza_Secundaria
    ## 68  Criollo blanco o pilchipahuak      Olotillo              NA
    ## 72   MaÌ_z blanco hueyi chipahuak      Olotillo              NA
    ## 75                 Olotillo crema      Olotillo              NA
    ## 79                 MaÌ_z olotillo      Olotillo              NA
    ## 123                        Chimbo      Olotillo              NA
    ## 159                      Olotillo      Olotillo              NA
    ##     AÌ.o._de_colecta                    Localidad    Municipio Estado.1
    ## 68              2009                Ixtlahuatempa     Huejutla  Hidalgo
    ## 72              2006                    Tenexco I    Atlapexco  Hidalgo
    ## 75              2009       Nuevo Vicente Guerrero   Villacorzo  Chiapas
    ## 79              2009              Ejido El Parral   Villacorzo  Chiapas
    ## 123             2009 Primera SecciÌ_n Medio Monte Tuxtla Chico  Chiapas
    ## 159             2009                  Medio Monte Tuxtla Chico  Chiapas
    ##      Longitud  Latitud Altitud     Ruizetal2008_grupo   Sanchezetal_grupo
    ## 68  -98.43611 21.13333     225 3A_muycalido990-1360mm MaduraciÌ_n tardÌ_a
    ## 72  -98.34639 21.06389     194 3A_muycalido990-1360mm MaduraciÌ_n tardÌ_a
    ## 75  -92.97972 16.03225     618 3A_muycalido990-1360mm MaduraciÌ_n tardÌ_a
    ## 79  -93.00194 16.36403     654 3A_muycalido990-1360mm MaduraciÌ_n tardÌ_a
    ## 123 -92.21693 14.90176     191 3A_muycalido990-1360mm MaduraciÌ_n tardÌ_a
    ## 159 -92.20450 14.87376     159 3A_muycalido990-1360mm MaduraciÌ_n tardÌ_a
    ##     Categ.Altitud            Rbiogeo              DivFloristic PeralesBiog
    ## 68            low    Golfo de Mexico Costa del Golfo de Mexico Planicies d
    ## 72            low    Golfo de Mexico Costa del Golfo de Mexico Planicies d
    ## 75            low Costa del Pacifico    Serranias Transismicas     Chiapas
    ## 79            low Costa del Pacifico            Costa Pacifica     Chiapas
    ## 123           low Costa del Pacifico                 Soconusco     Chiapas
    ## 159           low Costa del Pacifico                 Soconusco     Chiapas

## Pregunta 10 Crear un nuevo data frame solo con las muestras de las razas Reventador, Jala y Ancho

    subrazas <- subset(maiz, Raza %in% c("Reventador", "Jala", "Ancho"))
    head(subrazas)
    
    ##     OrderColecta NSiembra      Origen       Raza            Estado Num_Colecta
    ## 16             1       27    2009-232      Ancho Estado de MÌ©xico    2009-232
    ## 30             2       47 Ig-2010-294       Jala           Nayarit      2003-2
    ## 37             1       30    2009-236      Ancho Estado de MÌ©xico    2009-236
    ## 43             2       46 Ig-2010-201       Jala           Nayarit     NAY-129
    ## 57             1       29    2009-235      Ancho Estado de MÌ©xico    2009-235
    ## 103            1      114     2007-09 Reventador            Sonora SON2007-009
    ##       Nombre_comun Raza_Primaria Raza_Secundaria AÌ.o._de_colecta
    ## 16           Ancho         Ancho              NA             2009
    ## 30         Criollo          Jala              NA             2003
    ## 37           Ancho         Ancho              NA             2009
    ## 43         Criollo          Jala              NA             1951
    ## 57  Ancho Amarillo         Ancho              NA             2009
    ## 103     Reventador    Reventador              NA             2007
    ##                         Localidad    Municipio          Estado.1   Longitud
    ## 16  Carretera San Juan Tepeoculco     Atlautla Estado de MÌ©xico  -98.77861
    ## 30                        Jomulco         Jala           Nayarit -104.42869
    ## 37           San AndrÌ©s Tlalamac     Atlautla Estado de MÌ©xico  -98.80944
    ## 43                           Jala         Jala           Nayarit -104.44056
    ## 57  Carretera San Juan Tepeoculco     Atlautla Estado de MÌ©xico  -98.77861
    ## 103                     La Isleta Ì\u0081lamos            Sonora -108.91194
    ##      Latitud Altitud              Ruizetal2008_grupo Sanchezetal_grupo
    ## 16  18.99750    2226                          4_Jala      Ocho Hileras
    ## 30  21.10167    1153 2A_semicalido500-870mm y 4_Jala      Ocho Hileras
    ## 37  18.96806    2073                          4_Jala      Ocho Hileras
    ## 43  21.10000    1060 2A_semicalido500-870mm y 4_Jala      Ocho Hileras
    ## 57  18.99750    2226                          4_Jala      Ocho Hileras
    ## 103 26.84361     204          2A_semicalido500-870mm         Chapalote
    ##     Categ.Altitud            Rbiogeo           DivFloristic PeralesBiog
    ## 16           high      Eje Volcanico Serranias Meridionales Mesa Centra
    ## 30            mid Costa del Pacifico Serranias Meridionales Sierras del
    ## 37           high      Eje Volcanico Serranias Meridionales Mesa Centra
    ## 43            mid Costa del Pacifico Serranias Meridionales Sierras del
    ## 57           high      Eje Volcanico Serranias Meridionales Mesa Centra
    ## 103           low Costa del Pacifico         Costa Pacifica Sierras del

## Pregunta 11

    write.csv(subrazas,"C:/Users/pame5/Documents/BioinfinvRepro/Unidad1/Sesion3/PracUni1Ses3/maices/meta/submat.csv", row.names = FALSE)
    getwd()
    
    ## [1] "C:/Users/pame5/OneDrive/Escritorio/BIOINFORMATICA"
