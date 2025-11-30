Туториал
========

Этот туториал покажет, как использовать библиотеку ORCA Descriptors как через Python, так и через командную строку.

Использование как Python библиотека
-------------------------------------

Базовое использование
~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~

Библиотека автоматически кеширует результаты расчетов. Если вы рассчитываете дескрипторы для той же молекулы с теми же параметрами, будет использован кешированный результат::

   # Первый расчет - запускает ORCA
   homo1 = orca.homo_energy(mol)  # Занимает время
   
   # Второй расчет - использует кеш
   homo2 = orca.homo_energy(mol)  # Мгновенно

Использование как консольная утилита
-------------------------------------

Библиотеку также можно использовать как консольную утилиту после установки.

Запуск бенчмарка
~~~~~~~~~~~~~~~~

Перед оценкой времени расчетов запустите бенчмарк::

   orca_descriptors run_benchmark --working_dir ./calculations

Оценка времени расчета
~~~~~~~~~~~~~~~~~~~~~~

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
------------------------

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~

- ``get_atom_charges(mol)`` - Атомные заряды Малликена
- ``get_min_h_charge(mol, method="ESP")`` - Минимальный заряд водорода

Геометрические дескрипторы
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``xy_shadow(mol)`` - Площадь проекции на плоскость XY (Å²)
- ``molecular_volume(mol)`` - Молекулярный объем (Å³)
- ``get_bond_lengths(mol, atom1, atom2)`` - Длины связей (Å)

Дескрипторы реакционной способности
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``meric(mol)`` - Минимальный индекс электрофильности для углерода (эВ)
- ``dipole_moment(mol)`` - Дипольный момент (Дебай)

Топологические дескрипторы
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``topological_distance(mol, atom1, atom2)`` - Сумма топологических расстояний
- ``num_rotatable_bonds(mol)`` - Количество вращающихся связей
- ``wiener_index(mol)`` - Индекс Винера

Физико-химические дескрипторы
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``m_log_p(mol)`` - Коэффициент распределения Моригучи (октанол/вода)
- ``polar_surface_area(mol)`` - Полярная площадь поверхности (Å²)
- ``solvent_accessible_surface_area(mol)`` - Доступная площадь поверхности растворителя (Å²)

Термодинамические дескрипторы
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``gibbs_free_energy(mol)`` - Энергия Гиббса (Хартри)
- ``entropy(mol)`` - Энтропия (Дж/(моль·К))
- ``enthalpy(mol)`` - Энтальпия (Хартри)

Дескрипторы автокорреляции
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``moran_autocorrelation(mol, lag, weight)`` - Автокорреляция Моран
- ``autocorrelation_hats(mol, lag, unweighted)`` - Автокорреляция HATS

