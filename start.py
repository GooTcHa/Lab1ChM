import copy
from math import sqrt


def create_matrix(a, b, e):
    #Заполняем матрицы по правилам, преведённым в условии
    for i in range(1, 16):
        e.append([])
        a.append([])
        b.append(3*sqrt(i))
        for j in range(1, 16):
            if i == j:
                e[i - 1].append(1)
                a[i - 1].append(5*i)
            else:
                e[i - 1].append(0)
                a[i - 1].append(-(i+sqrt(j)))
    return a, b

#all
def get_left_diag(a, b, e):
    frac: float
    # print(*a, sep='\n')
    for k in range(0, 14):
        #Каждый цикл с i образует нулевой столбец под [k][k] элементом матрицы
        for i in range(k, 14):
            #Находим делитель
            frac = a[i + 1][k] / a[k][k]
            b[i + 1] -= frac * b[k]
            for j in range(0, 15):
                if j == k:
                    a[i + 1][j] = 0
                else:
                    #Находим разницу соответсвующих координат в строке k и i + k
                    a[i + 1][j] -= frac * a[k][j]
                #ПО методу Гаусса высчитываем обратную матрицы(на этом этапе приводя её к нижнедиагональной)
                e[i + 1][j] -= frac * e[k][j]
    return a, b, e


def get_right_diag(a, b, e):
    frac: float
    for k in range(0, 14):
        #Окончание этого цикла означает, что над a[14 - k][14 - k] нули
        for i in range(k, 14):
            #Находим делитель
            frac = a[13 - i][14 - k] / a[14 - k][14 - k]
            b[13 - i] -= frac * b[14 - k]
            for j in range(0, 15):
                if j == 14 - k:
                    a[13 - i][j] = 0
                else:
                    #Отнимаем соответствующие элементы с строках
                    a[13 - i][j] -= frac * a[14 - k][j]
                e[13 - i][j] -= frac * e[14 - k][j]
    return a, b


def get_answer(a, b, e):
    ans = []
    for i in range(15):
        #Находим ответ деля элемент вектора b на соответствующий элемент диагонали A
        ans.append(b[i]/a[i][i])
        #Получаем обратную матрицу
        for j in range(15):
            e[i][j] /= a[i][i]
    return ans


def get_residual(a, x, b):
    residual = []
    #Подставляем ответ в Ax=b и получаем невязку
    for i in range(14):
        for j in range(14):
            b[i] -= a[i][j] * x[j]
        residual.append(b[i])
    return residual


def multiply_matrix(a, b, e):
    # result = []
    s = 0.0
    m = -100000000000.0
    #Считаем ||E - AA^-1||_1
    for i in range(15):
        # result.append([])
        for j in range(15):
            s = 0
            for k in range(15):
                s += a[i][k] * b[k][j]
            # result[i].append(s)
            e[i][j] -= s
        m = max(sum(e[i]), m)
    return m


if __name__ == "__main__":
    a = []
    b = []
    e = []
    #Создаём матрицы
    create_matrix(a, b, e)
    #Выводим полученные матрицы
    print('Matrix A: ', *a, 'Matrix B: ', *b, 'E: ', *e, sep='\n')
    #Копируем матрицы, для последующей проверки ответов
    aCh = copy.deepcopy(a)
    bCh = copy.deepcopy(b)
    eCh = copy.deepcopy(e)
    #Приводим матрицу А к нижнедиагонильной форме
    get_left_diag(a, b, e)
    print('Diag A: ', *a, 'Matrix B: ', *b, 'A^-1: ', *e, sep='\n')
    #Приводим матрицу А к диагональной форме
    get_right_diag(a, b, e)
    #Приводим матрицу А к единичной форме и получаем ответ
    print('Diag A: ', *a, 'Matrix B: ', *b, 'A^-1: ', *e, sep='\n')
    ans = get_answer(a, b, e)
    print('A^-1', *e, sep='\n')
    #Выводим ответ
    print('X = ', ans, sep='\n')
    #Выводим норму невязки
    print('Норма невязки: ', max(get_residual(aCh, ans, bCh)), sep='\n')
    #Вторая проверка на отклонение ответа
    print('||E - A*A^-1||_1 = ', multiply_matrix(a, e, eCh), sep='\n')
