# Tree_generator

Описание:
Генератор деревьев, написанный на языке Python, позволяет по заданным начальным данным получить случайное коалесцентное дерево. Генератор позволяет выводить историю миграций и коалесценций в бинарный и в текстовый файлы в спецефическом формате (подробнее см. ниже), а также визуализировать сгенерированные коалесцентные деревья.

Запуск:
При запуске через консоль генератор принимает следующие дополнительные аргументы:
1) -txt. При установлении -txt=True сгенерированная история коалесценций и миграций будет выведерна в текстовый файл с названием "tree_history" в удобном для чтения формате. Пример: 
 - Migration of 0 from 1 to 0 at 1499 dt;
    Данная страка говорит о том, что во время t=1799 произошла миграция образца с Id = 0 из популяции №1 в популяцию № 0.
 - Coalescence in 2 population between 1 and 6; at 3534 dt;
    Данная строка говорит о том, что во время t=3534 в популяции №2 произошла коалесценция между образцами с Id 1 и 6.

2) -ac. При установлении аргумента -ac=n, где n - целое положительное значение будет сгенерированно n различных коалесцентных деревьев с одинаковыми начальными данными и посчитанно средний TMRCA (time to the most recent common ancestor) для этих деревьев. Среднее значение TMRCA выводится в консоль.

3) -d. При становлении -d=True сгенерированное дерево будет визуализироанно, картинка откроется в новом окне.
