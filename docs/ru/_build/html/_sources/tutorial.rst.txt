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

