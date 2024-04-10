--- 
title: "Diplomová práce"
author: "Tomáš Pompa"
date: "10. 04. 2024"
site: bookdown::bookdown_site
---

# SVM pro funkcionální data {-}

Cílem bude aplikovat poznatky o metodě podpůrných vektorů (SVM) pro mnohorozměrná data na data funkcionálního typu, tedy nekonečně-rozměrné objekty.
K tomu využijeme jednak převod (redukci) objektů z nekonečné dimenze na objekty konečné dimenze a následným využitím známých postupů a také modifikaci SVM přímo pro funkcionální data, k čemuž využijeme poznatky o Hilbertových prostorech a skalárním součinu.

Dalším cílem bude porovnání jednotlivých metod pro klasifikaci funkcionálních dat na reálných a simulovaných datech. Bylo by dobré vymyslet nějakou zajímavou simulační studii, která bude demonstrovat různá chování uvažovaných metod.

Mezi uvažované klasifikační metody patří:

- $K$ nejbližších sousedů (KNN),

- logistická regrese (jak obyčejná (LR) tak její funkcionální modifikace (LR_fda)),

- lineární (LDA) a kvadratická (QDA) diskriminační analýza,

- rozhodovací stromy (DT),

- náhodné lesy (RF) a 

- Support Vector Machines.

Postupně jednotlivé metody projdeme, nejprve na simulovaných datech, a následně budeme konstruovat metodu podpůrných vektorů pro funkcionální data (SVM_fda).

Základním balíčkem v `R` pro práci s funkcionálními objekty je `fda`. Dalšími užitečnými balíčky budou `MASS`, `e1071`, `fda.usc`, `refund` a další.
