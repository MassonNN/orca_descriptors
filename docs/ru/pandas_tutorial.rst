Туториал по интеграции с Pandas
============================

Этот туториал демонстрирует, как использовать библиотеку ORCA Descriptors вместе с pandas для эффективной пакетной обработки молекулярных дескрипторов.

Установка
------------

Для использования интеграции с pandas установите pandas как опциональную зависимость::

   pip install 'orca-descriptors[pandas]'

Или установите pandas отдельно::

   pip install pandas

Базовое использование
-----------

Метод ``calculate_descriptors`` обеспечивает бесшовную интеграцию с pandas DataFrame и Series.

Расчет всех дескрипторов
~~~~~~~~~~~~~~~~~~~~~~~~~~~

По умолчанию ``calculate_descriptors`` вычисляет все доступные дескрипторы для каждой молекулы::

   from orca_descriptors import Orca
   import pandas as pd

   # Инициализация калькулятора ORCA
   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )

   # Создание DataFrame со SMILES строками
   df = pd.DataFrame({
       'smiles': ['C1=CC=CC=C1', 'CCO', 'CC(=O)C'],
       'name': ['Бензол', 'Этанол', 'Ацетон']
   })

   # Расчет всех дескрипторов
   df_with_descriptors = orca.calculate_descriptors(smiles_column=df['smiles'])

   # DataFrame теперь содержит все колонки дескрипторов
   print(df_with_descriptors.columns)
   # Вывод: Index(['smiles', 'name', 'homo_energy', 'lumo_energy', 'gap_energy', ...])

   # Доступ к значениям дескрипторов
   print(df_with_descriptors[['name', 'homo_energy', 'lumo_energy', 'gap_energy']])

Выбор конкретных дескрипторов
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Вы можете указать, какие дескрипторы вычислять, используя параметр ``descriptors``::

   # Расчет только конкретных дескрипторов
   selected_descriptors = [
       'homo_energy',
       'lumo_energy',
       'gap_energy',
       'dipole_moment',
       'molecular_volume'
   ]

   df_subset = orca.calculate_descriptors(
       smiles_column=df['smiles'],
       descriptors=selected_descriptors
   )

   # Вычисляются только выбранные дескрипторы
   print(df_subset.columns)
   # Вывод: Index(['smiles', 'name', 'homo_energy', 'lumo_energy', 'gap_energy', 'dipole_moment', 'molecular_volume'])

Работа с Series
-------------------

Вы также можете передать pandas Series напрямую::

   # Создание Series со SMILES
   smiles_series = pd.Series(['C1=CC=CC=C1', 'CCO', 'CC(=O)C'])

   # Расчет дескрипторов
   df_result = orca.calculate_descriptors(smiles_column=smiles_series)

   # Результат - DataFrame с колонками 'smiles' и всеми дескрипторами
   print(df_result.head())

Отслеживание прогресса
-----------------

Метод по умолчанию предоставляет информацию о прогрессе::

   # Обработка множества молекул с выводом прогресса
   large_df = pd.DataFrame({
       'smiles': ['C1=CC=CC=C1'] * 10  # 10 молекул
   })

   df_result = orca.calculate_descriptors(
       smiles_column=large_df['smiles'],
       progress=True  # По умолчанию: показывает прогресс
   )

   # Вывод:
   # INFO - Processing molecule 1/10: C1=CC=CC=C1
   # INFO - Processing molecule 2/10: C1=CC=CC=C1
   # ...

Для отключения вывода прогресса::

   df_result = orca.calculate_descriptors(
       smiles_column=df['smiles'],
       progress=False
   )

Пример: Подготовка QSAR датасета
----------------------------------

Полный пример подготовки QSAR датасета::

   import pandas as pd
   from orca_descriptors import Orca

   # Инициализация калькулятора
   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )

   # Загрузка вашего молекулярного датасета
   df = pd.read_csv('molecules.csv')  # Содержит колонку 'smiles'

   # Расчет всех дескрипторов
   df_descriptors = orca.calculate_descriptors(smiles_column=df['smiles'])

   # Фильтрация молекул на основе дескрипторов
   active_molecules = df_descriptors[df_descriptors['gap_energy'] < 3.0]

   # Сохранение результатов
   df_descriptors.to_csv('molecules_with_descriptors.csv', index=False)

   # Статистический анализ
   print(df_descriptors[['homo_energy', 'lumo_energy', 'gap_energy']].describe())

Пример: Выбор пользовательских дескрипторов
-------------------------------------

Вычисление только энергетических дескрипторов для большого датасета::

   energy_descriptors = [
       'homo_energy',
       'lumo_energy',
       'gap_energy',
       'total_energy',
       'gibbs_free_energy',
       'entropy',
       'enthalpy'
   ]

   df_energy = orca.calculate_descriptors(
       smiles_column=df['smiles'],
       descriptors=energy_descriptors
   )

   # Анализ энергетических трендов
   print(df_energy.groupby('smiles')[energy_descriptors].mean())

Обработка ошибок
--------------

Метод корректно обрабатывает ошибки. Если расчет дескриптора для молекулы не удался, значение устанавливается в ``None``, и обработка продолжается::

   # Некоторые молекулы могут не пройти расчет дескрипторов
   df = pd.DataFrame({
       'smiles': ['C1=CC=CC=C1', 'INVALID_SMILES', 'CCO']
   })

   df_result = orca.calculate_descriptors(smiles_column=df['smiles'])

   # Неудачные расчеты помечены как None
   print(df_result[df_result['homo_energy'].isna()])

Доступные дескрипторы
---------------------

Следующие дескрипторы доступны для расчета:

**Энергетические дескрипторы:**
- ``homo_energy`` - Энергия HOMO (эВ)
- ``lumo_energy`` - Энергия LUMO (эВ)
- ``gap_energy`` - Зазор HOMO-LUMO (эВ)
- ``total_energy`` - Полная энергия (Хартри)
- ``gibbs_free_energy`` - Энергия Гиббса (Хартри)
- ``entropy`` - Энтропия (Дж/(моль·К))
- ``enthalpy`` - Энтальпия (Хартри)

**DFT дескрипторы:**
- ``ch_potential`` - Химический потенциал (эВ)
- ``electronegativity`` - Электроотрицательность (эВ)
- ``abs_hardness`` - Абсолютная жесткость (эВ)
- ``abs_softness`` - Абсолютная мягкость (1/эВ)
- ``frontier_electron_density`` - Максимальная плотность электронов на граничных орбиталях

**Дескрипторы полярной поверхности:**
- ``dipole_moment`` - Дипольный момент (Дебай)
- ``polar_surface_area`` - Полярная площадь поверхности (Å²)
- ``get_min_h_charge`` - Минимальный заряд водорода

**Термодинамические дескрипторы:**
- ``molecular_volume`` - Молекулярный объем (Å³)

**Многомерные дескрипторы:**
- ``num_rotatable_bonds`` - Количество вращающихся связей
- ``wiener_index`` - Индекс Винера
- ``solvent_accessible_surface_area`` - SASA (Å²)

**Геометрические дескрипторы:**
- ``xy_shadow`` - Площадь проекции на плоскость XY (Å²)

**Дескрипторы реактивности:**
- ``meric`` - Минимальный индекс электрофильности для углерода (эВ)

**Физико-химические дескрипторы:**
- ``m_log_p`` - Моригути Log P (коэффициент распределения октанол/вода)

**Дескрипторы автокорреляции:**
- ``moran_autocorrelation`` - Автокорреляция Моран (lag=2, weight='vdw_volume')
- ``autocorrelation_hats`` - Автокорреляция HATS (lag=4, unweighted=True)

Советы и лучшие практики
------------------------

1. **Кэширование**: Библиотека автоматически кэширует результаты расчетов. Если вы пересчитываете дескрипторы для тех же молекул, будут использованы кэшированные результаты.

2. **Пакетная обработка**: Для больших датасетов обрабатывайте молекулы пакетами для управления памятью и временем вычислений.

3. **Выбор дескрипторов**: Используйте параметр ``descriptors`` для расчета только нужных дескрипторов, что может значительно сократить время вычислений.

4. **Мониторинг прогресса**: Оставляйте ``progress=True`` (по умолчанию) для отслеживания длительных расчетов.

5. **Валидация данных**: Проверяйте значения ``None`` в результирующем DataFrame для выявления молекул, для которых не удался расчет дескрипторов.

Пример: Пакетная обработка
~~~~~~~~~~~~~~~~~~~~~~~~~~

   # Обработка пакетами
   batch_size = 10
   all_results = []

   for i in range(0, len(df), batch_size):
       batch = df.iloc[i:i+batch_size]
       batch_result = orca.calculate_descriptors(
           smiles_column=batch['smiles'],
           descriptors=['homo_energy', 'lumo_energy', 'gap_energy']
       )
       all_results.append(batch_result)

   # Объединение результатов
   final_df = pd.concat(all_results, ignore_index=True)

