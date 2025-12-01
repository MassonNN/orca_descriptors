Туториал
========

Этот туториал покажет, как использовать библиотеку ORCA Descriptors как через Python, так и через командную строку.

Использование как Python библиотека
--------------------------

Базовое использование
~~~~~~~~~~~

Сначала импортируйте необходимые классы::

   from orca_descriptors import Orca
   from rdkit.Chem import MolFromSmiles, AddHs

Инициализируйте калькулятор ORCA с нужными настройками::

   orca = Orca(
       script_path="orca",
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       dispersion_correction="D3BJ",
       solvation_model="COSMO(Water)",
       n_processors=8,
       pre_optimize=True,  # Включить предоптимизацию геометрии с MMFF94
   )

Создайте молекулу из SMILES строки::

   mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))  # Бензол

Рассчитайте дескрипторы::

   # Энергетические дескрипторы
   homo = orca.homo_energy(mol)
   lumo = orca.lumo_energy(mol)
   gap = orca.gap_energy(mol)
   
   # DFT дескрипторы
   mu = orca.ch_potential(mol)
   chi = orca.electronegativity(mol)
   eta = orca.abs_hardness(mol)
   
   # Термодинамические дескрипторы
   energy = orca.total_energy(mol)
   gibbs = orca.gibbs_free_energy(mol)
   
   # Дескрипторы молекулярных орбиталей
   homo_minus_1 = orca.mo_energy(mol, index=-2)  # Энергия HOMO-1
   
   # Дескрипторы зарядов
   min_h_charge = orca.get_min_h_charge(mol, method="ESP")  # Минимальный заряд водорода
   
   # Геометрические дескрипторы
   xy_area = orca.xy_shadow(mol)  # Площадь проекции на плоскость XY
   
   # Дескрипторы реакционной способности
   meric = orca.meric(mol)  # Индекс электрофильности для углерода
   
   # Топологические дескрипторы
   t_oo = orca.topological_distance(mol, 'O', 'O')  # Сумма расстояний O-O
   nrot = orca.num_rotatable_bonds(mol)  # Количество вращающихся связей
   wiener = orca.wiener_index(mol)  # Индекс Винера
   
   # Физико-химические дескрипторы
   logp = orca.m_log_p(mol)  # Коэффициент распределения октанол/вода
   sasa = orca.solvent_accessible_surface_area(mol)  # SASA
   
   # Дескрипторы автокорреляции
   mats2v = orca.moran_autocorrelation(mol, lag=2, weight='vdw_volume')
   hats4u = orca.autocorrelation_hats(mol, lag=4, unweighted=True)

Кеширование
~~~~~~~

Библиотека автоматически кеширует результаты расчетов. Если вы рассчитываете дескрипторы для той же молекулы с теми же параметрами, будет использован кешированный результат::

   # Первый расчет - запускает ORCA
   homo1 = orca.homo_energy(mol)  # Занимает время
   
   # Второй расчет - использует кеш
   homo2 = orca.homo_energy(mol)  # Мгновенно

Выбор функционалов и методов
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Библиотека поддерживает как DFT методы, так и полуэмпирические методы:

DFT методы
^^^^^^^^^^

Для расчетов с высокой точностью используйте DFT функционалы с базисными наборами::

   # Высокоточный DFT расчет
   orca_dft = Orca(
       functional="PBE0",           # Гибридный функционал
       basis_set="def2-TZVP",      # Тройной дзета базисный набор
       method_type="Opt",           # Оптимизация геометрии
       dispersion_correction="D3BJ", # Дисперсионная коррекция
       n_processors=8,
   )

Распространенные DFT функционалы:
- ``PBE0`` - Гибридный GGA функционал, хороший баланс точности и скорости
- ``B3LYP`` - Популярный гибридный функционал
- ``M06-2X`` - Meta-GGA функционал, хорош для термохимии
- ``ωB97X-D`` - Range-separated гибрид с дисперсией

Распространенные базисные наборы:
- ``def2-SVP`` - Маленький, быстрый (по умолчанию)
- ``def2-TZVP`` - Тройной дзета, более точный
- ``def2-QZVP`` - Четверной дзета, очень точный, но медленный

Полуэмпирические методы
^^^^^^^^^^^^^^^^^^^^^^^^

Для быстрых расчетов больших молекул используйте полуэмпирические методы::

   # Быстрый полуэмпирический расчет
   orca_semi = Orca(
       functional="AM1",            # Полуэмпирический метод
       method_type="SP",            # Одноточечный расчет (без оптимизации)
       n_processors=1,
       pre_optimize=True,           # Предоптимизировать с MMFF94
   )

Поддерживаемые полуэмпирические методы:
- ``AM1`` - Austin Model 1, хорош для органических молекул
- ``PM3`` - Parametric Method 3, улучшенная версия AM1
- ``PM6`` - Parametric Method 6, лучше для переходных металлов
- ``PM7`` - Parametric Method 7, улучшенная точность
- ``RM1`` - Recife Model 1, оптимизирован для органических соединений

Примечание: Для полуэмпирических методов параметры ``basis_set`` и ``dispersion_correction`` автоматически игнорируются.

Предоптимизация геометрии
~~~~~~~~~~~~~~~~~~~~~~~~~~

По умолчанию библиотека выполняет предоптимизацию геометрии с использованием силового поля MMFF94 из RDKit перед отправкой молекулы в ORCA. Это может значительно ускорить расчеты ORCA, особенно для оптимизации геометрии::

   # С предоптимизацией (по умолчанию)
   orca = Orca(
       functional="PBE0",
       method_type="Opt",
       pre_optimize=True,  # По умолчанию: True
   )
   
   # Без предоптимизации
   orca = Orca(
       functional="PBE0",
       method_type="Opt",
       pre_optimize=False,
   )

Преимущества предоптимизации:
- Быстрее сходимость ORCA (требуется меньше шагов оптимизации)
- Более стабильные расчеты (лучшая стартовая геометрия)
- Снижение вычислительных затрат

Предоптимизация использует MMFF94, который требует явных атомов водорода. Библиотека автоматически добавляет водород, если необходимо, и генерирует 3D координаты, если они отсутствуют.

Типы расчетов
~~~~~~~~~~~~~

Выберите подходящий тип расчета в зависимости от ваших потребностей::

   # Одноточечный расчет энергии (самый быстрый)
   orca_sp = Orca(
       functional="PBE0",
       method_type="SP",  # Одноточечный расчет
   )
   
   # Оптимизация геометрии (медленнее, но дает оптимизированную геометрию)
   orca_opt = Orca(
       functional="PBE0",
       method_type="Opt",  # Оптимизация
   )

Примечание: Некоторые дескрипторы (например, ``molecular_volume``, ``polar_surface_area``, ``solvent_accessible_surface_area``) требуют оптимизированной геометрии и могут работать некорректно с ``method_type="SP"``.

Использование как консольная утилита
-----------------------------

Библиотеку также можно использовать как консольную утилиту после установки.

Запуск бенчмарка
~~~~~~~~~~~~~

Перед оценкой времени расчетов запустите бенчмарк::

   orca_descriptors run_benchmark --working_dir ./calculations

Оценка времени расчета
~~~~~~~~~~~~~~~~~~~~~~~~~

Оцените, сколько времени займет расчет::

   orca_descriptors approximate_time --molecule CCO --method_type Opt

Все параметры ORCA доступны как аргументы CLI. Например::

   orca_descriptors approximate_time \\
       --molecule CCO \\
       --functional PBE0 \\
       --basis_set def2-TZVP \\
       --n_processors 4 \\
       --method_type Opt

Пример рабочего процесса
----------------

Вот полный пример расчета дескрипторов для нескольких молекул::

   from orca_descriptors import Orca
   from rdkit.Chem import MolFromSmiles, AddHs
   
   # Инициализация калькулятора
   orca = Orca(
       working_dir="./calculations",
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )
   
   # Список молекул для обработки
   smiles_list = [
       "C1=CC=CC=C1",  # Бензол
       "CCO",           # Этанол
       "CC(=O)C",       # Ацетон
   ]
   
   results = []
   for smiles in smiles_list:
       mol = AddHs(MolFromSmiles(smiles))
       results.append({
           "smiles": smiles,
           "homo": orca.homo_energy(mol),
           "lumo": orca.lumo_energy(mol),
           "gap": orca.gap_energy(mol),
           "dipole": orca.dipole_moment(mol),
       })
   
   # Обработка результатов
   for r in results:
       print(f"{r['smiles']}: Gap = {r['gap']:.2f} eV")

Доступные дескрипторы
---------------------

Библиотека предоставляет комплексный набор дескрипторов для QSAR анализа:

Энергетические дескрипторы
~~~~~~~~~~~~~~~~~~

- ``homo_energy(mol)`` - Энергия HOMO (эВ)
- ``lumo_energy(mol)`` - Энергия LUMO (эВ)
- ``gap_energy(mol)`` - Разрыв HOMO-LUMO (эВ)
- ``mo_energy(mol, index)`` - Энергия молекулярной орбитали по индексу (эВ)
- ``total_energy(mol)`` - Полная энергия (Хартри)

DFT дескрипторы
~~~~~~~~~~~~~~~

- ``ch_potential(mol)`` - Химический потенциал (эВ)
- ``electronegativity(mol)`` - Электроотрицательность (эВ)
- ``abs_hardness(mol)`` - Абсолютная жесткость (эВ)
- ``abs_softness(mol)`` - Абсолютная мягкость (1/эВ)
- ``frontier_electron_density(mol)`` - Плотность фронтирных электронов

Дескрипторы зарядов
~~~~~~~~~~~~~~~~~~

- ``get_atom_charges(mol)`` - Атомные заряды Малликена
- ``get_min_h_charge(mol, method="ESP")`` - Минимальный заряд водорода

Геометрические дескрипторы
~~~~~~~~~~~~~~~~~~~~~

- ``xy_shadow(mol)`` - Площадь проекции на плоскость XY (Å²)
- ``molecular_volume(mol)`` - Молекулярный объем (Å³)
- ``get_bond_lengths(mol, atom1, atom2)`` - Длины связей (Å)

Дескрипторы реакционной способности
~~~~~~~~~~~~~~~~~~~~~~

- ``meric(mol)`` - Минимальный индекс электрофильности для углерода (эВ)
- ``dipole_moment(mol)`` - Дипольный момент (Дебай)

Топологические дескрипторы
~~~~~~~~~~~~~~~~~~~~~~~

- ``topological_distance(mol, atom1, atom2)`` - Сумма топологических расстояний
- ``num_rotatable_bonds(mol)`` - Количество вращающихся связей
- ``wiener_index(mol)`` - Индекс Винера

Физико-химические дескрипторы
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``m_log_p(mol)`` - Коэффициент распределения Моригучи (октанол/вода)
- ``polar_surface_area(mol)`` - Полярная площадь поверхности (Å²)
- ``solvent_accessible_surface_area(mol)`` - Доступная площадь поверхности растворителя (Å²)

Термодинамические дескрипторы
~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``gibbs_free_energy(mol)`` - Энергия Гиббса (Хартри)
- ``entropy(mol)`` - Энтропия (Дж/(моль·К))
- ``enthalpy(mol)`` - Энтальпия (Хартри)

Дескрипторы автокорреляции
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``moran_autocorrelation(mol, lag, weight)`` - Автокорреляция Моран
- ``autocorrelation_hats(mol, lag, unweighted)`` - Автокорреляция HATS

