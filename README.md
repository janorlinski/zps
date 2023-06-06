# zps

## Jak zrobic listę wszystkich plikow txt?

```
ls 100umRun11*.txt > lista.list'
```

To będzie lista wszystkich widm (wszystkie subruny, wszystkie detektory) dla Runu 11. Można przygotować też listę dla wszystkich subrunów dla wybranego detektora:

```
ls 100umRun11no*Ge01.txt > lista.list
```

Oczywiście nazwa listy może być dowolna, ale koniecznie musi kończyć się na ".list".

2. jak uruchomić program do przeprowadzania kalibracji?

Na neutronxie odpalić środowisko root6 wpisując

```
root6
```

Odpowiednią listę plików .txt należy wpisać w głównej funkcji (na samym dole) makra convertTextToHisto.C. Program najlepiej uruchomić z flagą '-b' i skompilować:

```
root -b convertTextToHisto.C++
```

3. Program zapisuje: wszystkie kalibracje w plikach .cal o nazwach odpowiadających plikom .txt, wszystkie widma w postaci histogramów oraz histogramy kontrolne z fitowania w pliku ROOTa o nazwie takiej samej jak lista (tylko kończącej się na .root). 


