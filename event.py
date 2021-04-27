import struct


class Event:                   # для event_type: 0 - коалесценция, 1 - миграция
    def __init__(self, event_type, time, population=None, l1=None, l2=None, migrant=None, source=None, dist=None):
        # тут нужно делать проверку, все ли необходимые переменные заданы в вызове
        self.__type = event_type
        self.__time = int(time)
        if event_type == 0:
            self.__population = population
            self.__lineage = [l1, l2]
        elif event_type == 1:
            self.__migrant = migrant
            self.__source = source
            self.__dist = dist

    @property
    def type(self):
        return self.__type

    @property
    def time(self):
        return self.__time

    @property
    def lineage(self):
        return self.__lineage

    @property
    def source(self):
        return self.__source

    @property
    def dist(self):
        return self.__dist

    @property
    def migrant(self):
        return self.__migrant

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # это метод для временного вывода в текстовый файл
    # чтобы можно было удобно проверять логику программы
    def __str__(self):
        """
        if self.__population == -1:
            line = "Coalescence in original population"
            line += " between " + str(self.__lineage[0]) + " and " + str(self.__lineage[1]) + ";"
            line += " at " + str(self.__time) + " dt" + "\n"
        """
        if self.__type == 0:
            line = "Coalescence in " + str(self.__population) + " population"
            line += " between " + str(self.__lineage[0]) + " and " + str(self.__lineage[1]) + ";"
            line += " at " + str(self.__time) + " dt" + "\n"
        elif self.__type == 1:
            line = "Migration of " + str(self.__migrant) + " from " + str(self.__source) + " to " + str(self.__dist)
            line += " at " + str(self.__time) + " dt" + "\n"
        return line

    def to_binary_form(self):
        if self.__type == 0:              # если собитие - коалесценция
            return struct.pack("<5i", 0, self.__population, self.__lineage[0], self.__lineage[1], self.__time)
        elif self.__type == 1:            # если собитие - миграция
            return struct.pack("<5i", 1, self.__source, self.__dist, self.__migrant, self.__time)
