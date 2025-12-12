Туториал по пакетной обработке
================================

Этот туториал демонстрирует использование библиотеки ORCA Descriptors для эффективной пакетной обработки молекулярных дескрипторов с поддержкой интеграции с pandas, параллелизации и продвинутого определения дескрипторов.

Установка
---------

Для использования пакетной обработки с интеграцией pandas установите pandas как опциональную зависимость::

   pip install 'orca-descriptors[pandas]'

Или установите pandas отдельно::

   pip install pandas

Базовое использование
---------------------

Класс ``ORCABatchProcessing`` предоставляет эффективную пакетную обработку молекулярных дескрипторов.

Создание процессора пакетной обработки
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Инициализируйте процессор пакетной обработки с вашей конфигурацией ORCA::

   from orca_descriptors import Orca, ORCABatchProcessing

   # Создание экземпляра Orca
   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )

   # Создание процессора пакетной обработки
   batch_processing = ORCABatchProcessing(
       orca=orca,
       working_dir=".",
   )

Вы также можете создать процессор пакетной обработки без существующего экземпляра Orca::

   batch_processing = ORCABatchProcessing(
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )

Расчет дескрипторов
~~~~~~~~~~~~~~~~~~~

Рассчитайте дескрипторы для списка SMILES строк::

   smiles_list = ["C1=CC=CC=C1", "CCO", "CC(=O)C"]
   
   # Расчет всех доступных дескрипторов
   result = batch_processing.calculate_descriptors(smiles_list)
   
   # Результат - DataFrame (если pandas доступен) или список словарей
   print(result)

Работа с Pandas
---------------

Процессор пакетной обработки легко интегрируется с pandas DataFrames и Series.

Ввод DataFrame
~~~~~~~~~~~~~~

Передайте DataFrame с колонкой 'smiles'::

   import pandas as pd
   
   df = pd.DataFrame({
       'smiles': ['C1=CC=CC=C1', 'CCO', 'CC(=O)C'],
       'name': ['Бензол', 'Этанол', 'Ацетон']
   })
   
   # Расчет дескрипторов - исходные колонки сохраняются
   df_result = batch_processing.calculate_descriptors(df['smiles'])
   
   # DataFrame содержит исходные колонки + колонки дескрипторов
   print(df_result.columns)
   # Вывод: Index(['smiles', 'name', 'homo_energy', 'lumo_energy', ...])

Ввод Series
~~~~~~~~~~~

Вы также можете передать pandas Series напрямую::

   smiles_series = pd.Series(['C1=CC=CC=C1', 'CCO', 'CC(=O)C'])
   
   df_result = batch_processing.calculate_descriptors(smiles_series)
   
   # Результат - DataFrame только с колонками дескрипторов
   print(df_result.head())

Ввод списка
~~~~~~~~~~~

Поддерживаются обычные списки Python::

   smiles_list = ['C1=CC=CC=C1', 'CCO', 'CC(=O)C']
   
   result = batch_processing.calculate_descriptors(smiles_list)
   
   # Результат - DataFrame (если pandas доступен) или список словарей

Определение дескрипторов с помощью XMolecule API
-------------------------------------------------

XMolecule API позволяет определять дескрипторы с их параметрами, используя специальную молекулу-заглушку. Это особенно полезно для дескрипторов, требующих параметров.

Создание X молекулы
~~~~~~~~~~~~~~~~~~~

Создайте X молекулу с помощью метода ``x_molecule()``::

   x = batch_processing.x_molecule()

Использование X молекулы для определения дескрипторов
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Вызовите методы дескрипторов на экземпляре Orca с X молекулой для определения дескрипторов с параметрами::

   descriptors = [
       orca.ch_potential(x),
       orca.electronegativity(x),
       orca.abs_hardness(x),
       orca.topological_distance(x, 'O', 'O'),  # Расстояние между атомами кислорода
       orca.mo_energy(x, -3),  # Энергия HOMO-2 (индекс -3)
   ]
   
   result = batch_processing.calculate_descriptors(
       smiles_list,
       descriptors=descriptors
   )
   
   # Колонки результата включают имена дескрипторов с параметрами:
   # 'ch_potential', 'electronegativity', 'abs_hardness', 
   # 'topological_distance_O_O', 'mo_energy_-3'

Полный пример
~~~~~~~~~~~~~

Вот полный пример использования XMolecule API::

   from orca_descriptors import Orca, ORCABatchProcessing
   import pandas as pd
   
   # Инициализация
   orca = Orca(functional="PBE0", basis_set="def2-SVP")
   batch_processing = ORCABatchProcessing(orca=orca)
   
   # Создание X молекулы
   x = batch_processing.x_molecule()
   
   # Определение дескрипторов с параметрами
   descriptors = [
       orca.homo_energy(x),
       orca.lumo_energy(x),
       orca.gap_energy(x),
       orca.mo_energy(x, -1),  # HOMO
       orca.mo_energy(x, -2),  # HOMO-1
       orca.topological_distance(x, 'C', 'C'),  # Расстояния C-C
       orca.topological_distance(x, 'O', 'O'),  # Расстояния O-O
   ]
   
   # Загрузка вашего набора данных
   df = pd.read_csv('molecules.csv')
   
   # Расчет дескрипторов
   df_result = batch_processing.calculate_descriptors(
       df['smiles'],
       descriptors=descriptors
   )
   
   # Сохранение результатов
   df_result.to_csv('molecules_with_descriptors.csv', index=False)

Выбор дескрипторов
------------------

Вы можете указать, какие дескрипторы вычислять двумя способами:

Способ 1: Использование имен дескрипторов
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Передайте список имен дескрипторов как строки::

   selected_descriptors = [
       'homo_energy',
       'lumo_energy',
       'gap_energy',
       'dipole_moment',
       'molecular_volume'
   ]
   
   result = batch_processing.calculate_descriptors(
       smiles_list,
       descriptors=selected_descriptors
   )

Способ 2: Использование XMolecule API (Рекомендуется)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Используйте XMolecule API для дескрипторов с параметрами::

   x = batch_processing.x_molecule()
   
   descriptors = [
       orca.homo_energy(x),
       orca.mo_energy(x, -3),  # С параметром
       orca.topological_distance(x, 'O', 'O'),  # С параметрами
   ]
   
   result = batch_processing.calculate_descriptors(
       smiles_list,
       descriptors=descriptors
   )

Параллелизация
--------------

Процессор пакетной обработки поддерживает несколько режимов параллелизации для эффективной обработки больших наборов данных.

Последовательная обработка
~~~~~~~~~~~~~~~~~~~~~~~~~~

Режим по умолчанию - обработка молекул по одной::

   batch_processing = ORCABatchProcessing(
       orca=orca,
       parallel_mode="sequential"
   )

Мультипроцессинг
~~~~~~~~~~~~~~~~

Используйте Python multiprocessing для параллельного запуска нескольких расчетов ORCA::

   batch_processing = ORCABatchProcessing(
       orca=orca,
       parallel_mode="multiprocessing",
       n_workers=4  # Количество параллельных воркеров
   )
   
   # Обработка молекул параллельно
   result = batch_processing.calculate_descriptors(smiles_list)

Режим multiprocessing автоматически корректирует оценки времени на основе количества воркеров и эффективности параллелизации.

MPI (mpirun)
~~~~~~~~~~~~

Используйте встроенную параллелизацию MPI ORCA::

   batch_processing = ORCABatchProcessing(
       orca=orca,
       parallel_mode="mpirun",
       use_mpirun=True,
       n_processors=4  # Процессоры на расчет ORCA
   )

Отслеживание прогресса
----------------------

Процессор пакетной обработки предоставляет подробную информацию о прогрессе.

Вывод прогресса
~~~~~~~~~~~~~~~

По умолчанию отображается информация о прогрессе::

   result = batch_processing.calculate_descriptors(
       smiles_list,
       progress=True  # По умолчанию
   )
   
   # Вывод:
   # INFO - Estimating calculation times...
   # INFO - Processing molecule 1/10 (remaining: 10, estimated time: ~5m)
   # INFO - Processing molecule 2/10 (remaining: 9, estimated time: ~4m 30s, avg: 30.5s/molecule)
   # INFO - Processing molecule 3/10 (remaining: 8, CACHED, avg: 0.1s/molecule)

Кешированные молекулы
~~~~~~~~~~~~~~~~~~~~~~

Когда молекула найдена в кеше, прогресс показывает "CACHED" вместо оценки времени::

   # Первый расчет
   result1 = batch_processing.calculate_descriptors(smiles_list)
   
   # Второй расчет (из кеша)
   result2 = batch_processing.calculate_descriptors(smiles_list)
   # Вывод: INFO - Processing molecule 1/10 (remaining: 10, CACHED)

Отключение прогресса
~~~~~~~~~~~~~~~~~~~~

Для отключения вывода прогресса::

   result = batch_processing.calculate_descriptors(
       smiles_list,
       progress=False
   )

Обработка ошибок
---------------

Процессор пакетной обработки корректно обрабатывает ошибки. Если расчет дескриптора для молекулы не удался, значение устанавливается в ``None`` и обработка продолжается::

   df = pd.DataFrame({
       'smiles': ['C1=CC=CC=C1', 'INVALID_SMILES', 'CCO']
   })
   
   result = batch_processing.calculate_descriptors(df['smiles'])
   
   # Неудачные расчеты помечаются как None
   print(result[result['homo_energy'].isna()])

Сообщения об ошибках логируются на разных уровнях:

- Краткое резюме ошибки: ``logging.INFO``
- Подробная информация об ошибке: ``logging.DEBUG``

Пример: Подготовка QSAR набора данных
---------------------------------------

Вот полный пример подготовки QSAR набора данных::

   import pandas as pd
   from orca_descriptors import Orca, ORCABatchProcessing
   
   # Инициализация калькулятора
   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )
   
   # Создание процессора пакетной обработки с мультипроцессингом
   batch_processing = ORCABatchProcessing(
       orca=orca,
       parallel_mode="multiprocessing",
       n_workers=4,
   )
   
   # Создание X молекулы для определения дескрипторов
   x = batch_processing.x_molecule()
   
   # Определение дескрипторов
   descriptors = [
       orca.homo_energy(x),
       orca.lumo_energy(x),
       orca.gap_energy(x),
       orca.ch_potential(x),
       orca.electronegativity(x),
       orca.abs_hardness(x),
       orca.dipole_moment(x),
       orca.molecular_volume(x),
   ]
   
   # Загрузка вашего набора данных молекул
   df = pd.read_csv('molecules.csv')  # Содержит колонку 'smiles'
   
   # Расчет дескрипторов
   df_descriptors = batch_processing.calculate_descriptors(
       smiles_column=df['smiles'],
       descriptors=descriptors
   )
   
   # Фильтрация молекул на основе дескрипторов
   active_molecules = df_descriptors[df_descriptors['gap_energy'] < 3.0]
   
   # Сохранение результатов
   df_descriptors.to_csv('molecules_with_descriptors.csv', index=False)
   
   # Статистический анализ
   print(df_descriptors[['homo_energy', 'lumo_energy', 'gap_energy']].describe())

Пример: Пакетная обработка с пользовательскими параметрами
------------------------------------------------------------

Обработка молекул пакетами с разными конфигурациями::

   # Обработка пакетами
   batch_size = 10
   all_results = []
   
   for i in range(0, len(df), batch_size):
       batch = df.iloc[i:i+batch_size]
       
       x = batch_processing.x_molecule()
       descriptors = [
           orca.homo_energy(x),
           orca.lumo_energy(x),
           orca.gap_energy(x)
       ]
       
       batch_result = batch_processing.calculate_descriptors(
           smiles_column=batch['smiles'],
           descriptors=descriptors
       )
       all_results.append(batch_result)
   
   # Объединение результатов
   final_df = pd.concat(all_results, ignore_index=True)

Удаленный кеш
-------------

Процессор пакетной обработки поддерживает удаленное кеширование через API, позволяя делиться результатами расчетов между разными машинами и пользователями. Доступен публичный сервер кеша по адресу ``https://api.orca-descriptors.massonnn.ru``.

Для использования удаленного кеша:

1. Зарегистрируйтесь на `orca-descriptors.massonnn.ru <https://orca-descriptors.massonnn.ru>`_ и выпустите API токен
2. Укажите API токен при инициализации экземпляра Orca::

   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       cache_api_token="ваш_api_токен",
   )
   
   batch_processing = ORCABatchProcessing(orca=orca)

Удаленный кеш работает прозрачно в пакетной обработке:
- Все молекулы проверяются на наличие кеша (локального и удаленного) перед началом расчетов
- Кешированные молекулы обрабатываются сразу и исключаются из дальнейших расчетов
- Новые расчеты автоматически загружаются на удаленный сервер (если у вас есть права на загрузку)
- Сообщения о прогрессе показывают, когда молекулы найдены в кеше

Режим только кеша
-----------------

Вы можете включить режим только кеша для использования только закешированных результатов без запуска расчетов ORCA::

   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       cache_api_token="ваш_api_токен",
       cache_only=True,  # Только использовать кеш, не запускать расчеты
   )
   
   batch_processing = ORCABatchProcessing(orca=orca)
   
   result = batch_processing.calculate_descriptors(smiles_list)

В режиме только кеша:
- Используются только закешированные результаты (локальные или удаленные)
- Молекулы, не найденные в кеше, возвращают ``None`` для дескрипторов
- Расчеты ORCA не выполняются
- Сообщения о прогрессе показывают предупреждения для молекул, не найденных в кеше
- Полезно для быстрого получения результатов из кеша без выполнения дорогостоящих расчетов

Советы и лучшие практики
-------------------------

1. **Кеширование**: Библиотека автоматически кеширует результаты расчетов. Пересчет дескрипторов для тех же молекул использует кешированные результаты, значительно ускоряя последующие запуски.

2. **Удаленный кеш**: Используйте удаленный кеш для обмена результатами между машинами и совместной работы с другими. Зарегистрируйтесь на `orca-descriptors.massonnn.ru <https://orca-descriptors.massonnn.ru>`_ для получения API токена.

3. **Мультипроцессинг**: Для больших наборов данных используйте ``parallel_mode="multiprocessing"`` с подходящим количеством воркеров (обычно равным количеству ядер CPU).

4. **Выбор дескрипторов**: Используйте XMolecule API для определения дескрипторов с параметрами, делая ваш код более читаемым и поддерживаемым.

5. **Мониторинг прогресса**: Оставляйте ``progress=True`` (по умолчанию) для мониторинга длительных расчетов и идентификации кешированных молекул.

6. **Валидация данных**: Проверяйте значения ``None`` в результирующем DataFrame для идентификации молекул, для которых расчет дескрипторов не удался или которые не найдены в кеше (при использовании режима только кеша).

7. **Управление памятью**: Для очень больших наборов данных обрабатывайте молекулы пакетами для управления использованием памяти.

Доступные дескрипторы
---------------------

См. документацию :doc:`descriptors` для полного списка доступных дескрипторов и их параметров.

