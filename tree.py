import numpy as np
from scipy.special import binom  # для биномиальных коэффициентов
import functions as fun


class Tree:

    def __init__(self, coef_file: str, migr_file: str):
        # Извлечение данных из файла и создание объекта
        # Данные в файле должны быть в таком порядке: N, k, n, t, m, q, Q
        N, k, n, T, q, Q = fun.read_coef(coef_file)
        self.__N = int(N[0])
        self.__number_of_populations = int(k[0])
        self.__number_of_samples = np.array(n, dtype=int)
        self.__T = int(T[0])
        self.__migration_probability = np.array(fun.read_migration(migr_file), ndmin=2)
        self.__coalescence_probability = np.array(q, dtype=float)
        self.__Q = float(Q[0])

    # для отладки
    def show(self):
        print(self.__number_of_populations, self.__T, '\n')
        print(self.__number_of_samples, '\n')
        print(self.__migration_probability, '\n\n', self.__coalescence_probability, '\n')

    @property
    def original_size(self):
        return self.__N

    @property
    def number_of_populations(self):
        return self.__number_of_populations

    @property
    def number_of_samples(self):
        return self.__number_of_samples

    @property
    def T(self):
        return self.__T

    @property
    def migration_probability(self):
        return self.__migration_probability

    @property
    def coalescence_probability(self):
        return self.__coalescence_probability

    @property
    def Q(self):
        return self.__Q

    ####################################################
    # расчет вероятности события при данных размерах популяций
    # принимает список текущих размеров всех популяций
    def event_rate(self, population_sizes, migration=1):  # migration=1 - с учётом миграции, migration=0 - без учёта
        tmp_sum = 0
        # "вероятности" коалесценции
        for i in range(self.__number_of_populations):
            if population_sizes[i] <= 1:
                continue
            tmp_sum += binom(population_sizes[i], 2) * self.__coalescence_probability[i]
        # "вероятности" миграции
        if migration == 1:  # если миграция учитывается
            for i in range(self.__number_of_populations):  # по каждой популяции (из которой миграция)
                if population_sizes[i] == 0:
                    continue
                for j in range(self.__number_of_populations):  # по каждой популяции (в которую миграция)
                    if i == j:
                        continue
                    tmp_sum += self.__migration_probability[i][j] * population_sizes[i]
        rate = tmp_sum / (2 * self.__N)

        return rate

    #####################################################
    # расчет вероятности коалесценции в данной популяции
    # принимает номер популяции и её тикущий размер
    def count_coalescence_rate(self, pop_index, pop_size):
        result = self.__coalescence_probability[pop_index] * binom(pop_size, 2)
        return result / (2 * self.__N)

    #####################################################
    # расчет вероятности миграции из source_pop в dist_pop
    # принимает номера source_pop и dist_pop
    def count_migration_rate(self, source_pop_size, source_pop, dist_pop):
        result = source_pop_size * self.__migration_probability[source_pop][dist_pop]
        return result / (2 * self.__N)

    #####################################################
    # расчет вероятности коалесценции в original_population
    # принимает размер original_population
    def coalescence_rate_in_orig_pop(self, pop_size):
        result = self.__Q * binom(pop_size, 2)
        return result / (2 * self.__N)
