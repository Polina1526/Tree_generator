import numpy as np
from tree import Tree
from event import Event
import random
from struct import pack, unpack
import matplotlib.pyplot as plt


def read_coef(file_name):
    """
    :param file_name:
    :return: list of lists of coefficients in a certain order
    """
    with open(file_name) as f:
        line = next(f).rstrip()  # header with version etc
        line = line.split(" ")
        coefficients = []
        for line in f:
            if line[0] == "#":
                continue
            line = line.rstrip()
            line = line.split(" ")
            coefficients.append([float(value) for value in line])
        return coefficients


def read_migration(file_name):
    """
    :param file_name:
    :return: migration rates matrix (in format: list of lists)
    """
    with open(file_name) as f:
        line = next(f).rstrip()  # header with version etc
        line = line.split(" ")
        migration_rates = []
        for line in f:
            if line[0] == "#":
                continue
            line = line.rstrip()
            line = line.split(" ")
            migration_rates.append([float(value) for value in line])
        for i in range(len(migration_rates)):
            migration_rates[i][i] = 0.0
        return(migration_rates)


def generate_tree(tree: Tree, txt_file: bool) -> int:
    """
    :param tree: the initial data to generate a random tree
    :param txt_file: determines whether the tree history should be output to the txt file
    :return: nothing to return
    """
    ################################################
    # двумерный 'массив'(список списков), каждая строка которого содержит поэлементно образцы каждой популяции
    # нумерация образцов сквозная (за это и отвечат tmp)
    current_population_sizes = np.zeros(tree.number_of_populations, dtype=int)
    tmp = 0
    samples = []
    for i in range(tree.number_of_populations):
        samples.append(list(range(tmp, tmp + tree.number_of_samples[i])))
        current_population_sizes[i] = tree.number_of_samples[i]
        tmp += tree.number_of_samples[i]

    # вывод в текстовый файл, если требуется
    if txt_file:
        out_file = open("tree_history", 'w', encoding="utf8")
        if out_file.closed:
            print("unable to open 'tree_history' file!")

    # для записи данных в бинарный файл
    out_file_bin = open('tree_history.bin', 'wb')
    if out_file_bin.closed:
        print("unable to open 'tree_history' file!")
    # в начале файла будет записано TMRCA, оставлю место для этого числа (unsigned short)
    out_file_bin.write(pack('<I', 0))

    timer = 0  # отсчет времени ведётся с единицы
    lineage_ind = []
    current_number_of_samples = np.sum(tree.number_of_samples)
    event_occurred = 0  # вспомогательный флаг

    if tree.T == -1:
        pop_merge_time = timer + 1
    else:
        pop_merge_time = tree.T

    # если популяций больше одной
    if tree.number_of_populations > 1:
        # пока не останется всего один образец и пока популяции не объединились
        while current_number_of_samples > 1 and (timer + 1) <= pop_merge_time:
            timer += 1  # после увеличения показывает время текущей итерации
            if tree.T == -1:
                pop_merge_time = timer + 1
            rand = random.random()
            if rand >= tree.event_rate(current_population_sizes):  # если никаких событий в этом интервале не происходит
                continue
            # проверка коалесценции
            for i in range(tree.number_of_populations):
                if rand >= tree.count_coalescence_rate(i, current_population_sizes[i]):
                    # коалесценции в этой популяции не происходит
                    rand -= tree.count_coalescence_rate(i, current_population_sizes[i])
                else:
                    # происходит коалесценция в i-ой популяции:
                    # !!!!!!!!!!нужна подходящая библиотечная либо своя функция для этого места
                    # выбор двух случайных ИНДЕКСОВ образцов из sample[i]
                    if current_population_sizes[i] > 2:  # если нужно делать случайный выбор двух из всех
                        lineage_ind.append(random.randint(0, current_population_sizes[i] - 1))
                        lineage_ind.append(random.randint(0, current_population_sizes[i] - 2))
                        if lineage_ind[1] >= lineage_ind[0]:
                            lineage_ind[1] += 1
                    else:  # если образцов осталось всего два
                        lineage_ind.append(0)
                        lineage_ind.append(1)
                    lineage_ind.sort(
                        key=lambda lin: samples[i][lin])  # сортировка не по индексам, а по значениям ID элементов
                    # создаётся событие
                    event = Event(0, timer, population=i, l1=samples[i][lineage_ind[0]], l2=samples[i][lineage_ind[1]])
                    ################################################
                    if txt_file:  # информация о событии выводится в txt файл, если необходимо
                        out_file.write(str(event))
                    out_file_bin.write(event.to_binary_form())  # информация о событии выводится в бинарный файл
                    ################################################
                    # теперь из двух потомков останется один, с меньшем значением ID он и будет считаться предком
                    # на место потомка с большим ID встанет самый последний образец из всей популяции
                    if lineage_ind[1] == current_population_sizes[i] - 1:
                        samples[i].pop()  # удаление последнего элемента из списка
                    else:
                        samples[i][lineage_ind[1]] = samples[i][current_population_sizes[i] - 1]
                        samples[i].pop()  # удаление последнего элемента из списка
                    current_population_sizes[i] -= 1
                    current_number_of_samples -= 1
                    lineage_ind.clear()
                    # поднятие флага о произошедшем событии
                    event_occurred = 1
                    break

            # проверка, произошла ли коалесценция, есди да, то миграции не будет (переход к следующей итерации)
            if event_occurred == 1:
                event_occurred = 0
                continue

            # иначе проверка миграции
            for i in range(tree.number_of_populations):  # по каждой исходной популяции
                for j in range(tree.number_of_populations):  # по каждой возможной для миграции популяции
                    if i == j:  # пропускаем миграцию из i-ой популяции в i-ую
                        continue
                    if rand >= tree.count_migration_rate(current_population_sizes[i], i, j):
                        # миграции из i в j не происходит
                        rand -= tree.count_migration_rate(current_population_sizes[i], i, j)
                    else:  # миграция из i в j
                        # выбираем ИНДЕКС мигрирующего одразца
                        migrant_ind = random.randint(0, current_population_sizes[i] - 1)
                        # создаётся событие
                        event = Event(1, timer, migrant=samples[i][migrant_ind], source=i, dist=j)
                        ################################################
                        if txt_file:  # информация о событии выводится в txt файл, если необходимо
                            out_file.write(str(event))
                        out_file_bin.write(event.to_binary_form())  # информация о событии выводится в бинарный файл
                        ################################################
                        samples[j].append(samples[i][migrant_ind])
                        current_population_sizes[j] += 1
                        # на место мигранта встает последний образец данной популяции
                        if migrant_ind == current_population_sizes[i] - 1:  # если мигрант был "последним"
                            samples[i].pop()  # удаление последнего элемента из списка
                        else:
                            samples[i][migrant_ind] = samples[i][current_population_sizes[i] - 1]
                            samples[i].pop()  # удаление последнего элемента из списка
                        current_population_sizes[i] -= 1
                        # поднятие флага о произошедшем событии
                        event_occurred = 1
                        break
                # проверка, произошла ли миграция, есди да, то выход из первого цикла for (переход к следующей итерации)
                if event_occurred == 1:
                    event_occurred = 0
                    break

    # после достижения времени T популяции объединяются
    if current_number_of_samples > 1:  # если предок ещё не найден
        original_samples = []  # все "оставшиеся" образцы помещаются в один список
        for i in range(tree.number_of_populations):
            original_samples.extend(samples[i])

        # так как популяция теперь одна, то миграции нет, только коалесценция
        while current_number_of_samples > 1:
            timer += 1
            if random.random() <= tree.coalescence_rate_in_orig_pop(current_number_of_samples):
                # если коалесценция происходит:
                # !!!!!!!!!нужна подходящая библиотечная либо своя функция для этого места
                # выбор двух случайных ИНДЕКСОВ образцов из original_samples
                if current_number_of_samples > 2:  # если нужно делать случайный выбор двух из всех
                    lineage_ind.append(random.randint(0, current_number_of_samples - 1))
                    lineage_ind.append(random.randint(0, current_number_of_samples - 2))
                    if lineage_ind[1] >= lineage_ind[0]:
                        lineage_ind[1] += 1
                else:  # если образцов осталось всего два
                    lineage_ind.append(0)
                    lineage_ind.append(1)
                #################################
                lineage_ind.sort(
                    key=lambda lin: original_samples[lin])  # сортировка не по индексам, а по значениям ID элементов
                # создаётся событие(population=-1 значит что популяция общая)
                event = Event(0, timer, population=-1, l1=original_samples[lineage_ind[0]],
                              l2=original_samples[lineage_ind[1]])
                ################################################
                if txt_file:  # информация о событии выводится в txt файл, если необходимо
                    out_file.write(str(event))
                out_file_bin.write(event.to_binary_form())  # информация о событии выводится в файл
                ################################################
                # теперь из двух потомков останется один, с меньшем значением ID он и будет считаться предком
                # на место потомка с большим ID встанет самый последний образец из всей популяции
                if lineage_ind[1] == current_number_of_samples - 1:
                    original_samples.pop()  # удаление последнего элемента из списка
                else:
                    original_samples[lineage_ind[1]] = original_samples[current_number_of_samples - 1]
                    original_samples.pop()  # удаление последнего элемента из списка
                current_number_of_samples -= 1
                lineage_ind.clear()

    # запись TMRCA в начало файла
    out_file_bin.seek(0, 0)
    out_file_bin.write(pack('<I', timer))

    if txt_file:  # txt файл закрывается, если был открыт
        out_file.close()
    out_file_bin.close()

    return timer


def draw_tree(tree):
    out_file_bin = open("tree_history.bin", "rb")
    TMRCA, = unpack('<I', out_file_bin.read(4))
    N = np.sum(tree.number_of_samples)

    event_size = 4 * 5
    flag_time = (unpack('<5i', out_file_bin.read(event_size)))[4]  # время первого события

    # если слияние популяций было и было оно до TNRCA
    if tree.T != (-1) and tree.T < TMRCA:
        flag = 0
    else:
        flag = 1
    x = []
    y = []
    branch_width = []
    color = []  # текущий цвет ветви, ветвь обозначает текущую популяцию
    pop_change_ID = []  # фиксирует на каких моментах (индексы в масиве) менялась популяция
    for i in range(N):
        x.append([])
        y.append([])
        branch_width.append([0])
        color.append([])
        pop_change_ID.append([])
    color_scheme = {-1: "black",
                    0: "green",
                    1: "red",
                    2: "blue",
                    3: "yellow",
                    4: "magenta",
                    5: "cyan",
                    6: "white"}

    x[0].append(N / 2)
    y[0].append(TMRCA + 40)
    branch_width[0] = N

    # нахождение количества событий
    out_file_bin.seek(0, 2)
    file_size = out_file_bin.tell()
    event_size = 4 * 5
    event_amount = int((file_size - 4) / (4 * 5))  # так как в начале записано TMRCA

    # пока не прочитали все события из файла
    for i in range(event_amount):
        out_file_bin.seek(file_size - event_size * (i + 1))
        event = unpack('<5i', out_file_bin.read(event_size))

        # настало ли время слияния популяций в одну
        # flag нужен, чтобы не заходить в этот if второй раз
        if event[-1] < tree.T and flag == 0:
            for line in x:
                if len(line) > 0:
                    line.append(line[-1])
            for line in y:
                if len(line) > 0:
                    line.append(tree.T)
            for i in range(len(pop_change_ID)):
                if len(x[i]) > 0:
                    pop_change_ID[i].append(len(x[i]))
            flag = 1

        # если коалесценция
        if event[0] == 0:
            ev, pop_id, offspr1, offspr2, time = event
            if flag == 1 and len(color[offspr1]) > 0 and color[offspr1][-1] == color_scheme[-1]:
                color[offspr1].append(color_scheme[pop_id])
            # отметка точки коалесценции
            for offspr in [offspr1, offspr2]:
                x[offspr].append(x[offspr1][-1])  # т.к. 2й линии ещё нет, обе начинаются с x[offspr1][-1]
                y[offspr].append(time)
            # устанавливается цвет популяции
            if i == 0:
                for offspr in [offspr1, offspr2]:
                    color[offspr].append(color_scheme[pop_id])
            else:
                color[offspr2].append(color_scheme[pop_id])

            # раздвоение (отрисовка горизонтальной линии)
            x[offspr1].append(x[offspr1][-1] - branch_width[offspr1] / 4)
            x[offspr2].append(x[offspr2][-1] + branch_width[offspr1] / 4)
            for offspr in [offspr1, offspr2]:
                y[offspr].append(y[offspr][-1])  # y координата остаётся неизменной

            # изменение ширины ветви
            new_branch_width = branch_width[offspr1] / 2
            for offspr in [offspr1, offspr2]:
                branch_width[offspr] = new_branch_width

        # если миграция (при миграции меняется только цвет)
        if event[0] == 1:
            ev, source, dist, migrant, time = event

            if flag == 1 and len(color[migrant]) > 0 and color[migrant][-1] == color_scheme[-1]:
                color[migrant].append(color_scheme[dist])
            # добавляется точка, в которой была миграция
            x[migrant].append(x[migrant][-1])
            y[migrant].append(time)
            # запомниаю, на каком индексе (i) массива x[i] происходит смена цвета
            pop_change_ID[migrant].append(len(x[migrant]))

            # меняется цвет с dist -> source, так как дерево рисуется сверху вниз
            color[migrant].append(color_scheme[source])

    # проверка на то, что после слияния популяций все линии встречались в цикле и цвет был установлен
    for i in range(len(color)):
        if color[i][-1] == color_scheme[-1]:
            # установление цвета на цвет изначально популяции
            tmp_id = i
            for pop in range(tree.number_of_populations):
                # если не в этой популяции
                if tmp_id >= tree.number_of_samples[pop]:
                    tmp_id -= tree.number_of_samples[pop]
                else:
                    color[i].append(color_scheme[pop])
                    break

    # нужно опустить все ветви вниз до time = 0
    for line in x:
        line.append(line[-1])
    for line in y:
        line.append(0)

    for sample in range(N):
        pop_change_ID[sample].append(len(x[sample]))

    plt.figure(figsize=(3, 4))

    for sample in range(N):
        # если текущая ветка не меняла цвета
        if not pop_change_ID[sample]:  # если список пустой
            plt.plot(x[sample], y[sample], color=color[sample][0])

        else:
            start = 0
            for i in range(len(pop_change_ID[sample])):
                end = pop_change_ID[sample][i]
                plt.plot(x[sample][start:end], y[sample][start:end], color=color[sample][i], linewidth=0.9)
                start = end - 1

    if tree.T != -1:
        plt.axhline(y=tree.T, linestyle='--', linewidth=1)
    line = plt.scatter([N / 2], [TMRCA], color="black", label="o", s=5)

    plt.xticks([])
    plt.ylabel("Time in 2N populations")
    if tree.T == -1 or tree.T > TMRCA:
        handles = []
        labels = []
        color_for_legend = []
    else:
        handles = [line]
        labels = ["pooled pop"]
        color_for_legend = [color_scheme[-1]]
    for i in range(tree.number_of_populations):
        handles.append(line)
        labels.append(str(i + 1) + " pop")
        color_for_legend.append(color_scheme[i])
    plt.legend(handles=handles, labels=labels, labelcolor=color_for_legend)
    plt.show()
