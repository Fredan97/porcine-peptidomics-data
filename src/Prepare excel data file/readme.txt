Först körs filen "FixInExcel.py" i en folder som innehåller alla exported from peaks filerna från en specifik dag (folders som sedan innehåller protein-peptides.csv). Detta skapar en datafil för den dagen.

Dessa filer får sedan döpas om för till "Data Day 1.xlsx", "Data Day 2.xlsx" osv. 

Jag har sedan kört filerna "Fix excel file.py", följt av "Fix excel file 2.py". Nu finns en gemensam datafil för dag 1-3. 

När jag inkorporerade rerunsen så trimmade jag först ner den stora datafilen samt datafilen från rerunsen manuellt så att den inte innehöll kopplingar till specifika proteiner, och tog samtidigt bort duplikat av peptider som tillhörde flera olika protein. (Dvs. jag tog bort kolumnen som innehöll proteinnamnet, och tog sedan bort rader som var duplikat på avseende av peptidsekvensen). Jag döpte de filerna så att de slutar med OnlyPep NoDups och körde sedan de i "Fix excel file3.py" Nu är även de blindade proven inkorporerade men koppling till specifika protein är borta, vilket var relevant för den nya artikeln men inte för den gamla.