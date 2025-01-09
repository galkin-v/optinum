import numpy as np

def power_method(A: np.ndarray, num_iterations: int = 1000, tol: float = 1e-6) -> tuple:
    """
    Степенной метод для нахождения наибольшего по модулю собственного значения и соответствующего собственного вектора.

    Этот метод применяет итерационный процесс для нахождения наибольшего по модулю собственного значения матрицы A.
    На каждом шаге матрица умножается на вектор, и результат нормализуется. Метод прекращается, когда изменение
    собственного значения на текущем шаге становится меньше заданного порога.

    Параметры:
    -----------
    A : np.ndarray
        Квадратная матрица (размер n x n), для которой нужно найти собственное значение и собственный вектор.

    num_iterations : int, по умолчанию 1000
        Максимальное количество итераций для выполнения степенного метода.

    tol : float, по умолчанию 1e-6
        Порог сходимости. Метод завершится, если изменение собственного значения на текущем шаге
        будет меньше этого значения.

    Возвращает:
    --------
    tuple
        - собственное значение (float)
        - собственный вектор (np.ndarray), соответствующий наибольшему собственному значению

    Пример:
    --------
    >>> A = np.array([[4, 1],
                      [2, 3]])
    >>> lambda_max, eigenvector = power_method(A)
    >>> print(lambda_max)
    5.0
    >>> print(eigenvector)
    [0.57735027 0.57735027]

    В этом примере:
    - Матрица A = [[4, 1], [2, 3]]
    - Функция находит наибольшее по модулю собственное значение (5) и соответствующий собственный вектор.
    """

    # Инициализация начального вектора (случайный вектор)
    n = A.shape[0]
    x = np.random.rand(n)

    # Нормируем начальный вектор
    x = x / np.linalg.norm(x)

    # Инициализация переменной для хранения предыдущего собственного значения
    prev_lambda = 0

    for i in range(num_iterations):
        # Умножаем матрицу A на текущий вектор x
        x_next = np.dot(A, x)

        # Нормируем новый вектор
        x_next = x_next / np.linalg.norm(x_next)

        # Вычисляем приближенное собственное значение (скалярное произведение)
        lambda_ = np.dot(x_next.T, np.dot(A, x_next))

        # Проверяем сходимость (если изменение собственного значения меньше порога, выходим)
        if np.abs(lambda_ - prev_lambda) < tol:
            break

        # Обновляем значения для следующей итерации
        x = x_next
        prev_lambda = lambda_

    return lambda_, x_next

if __name__ == "__main__":
    A = np.array([[4, 1],
                  [2, 3]])
    lambda_max, eigenvector = power_method(A)
    print("Наибольшее собственное значение:", lambda_max)
    print("Соответствующий собственный вектор:", eigenvector)