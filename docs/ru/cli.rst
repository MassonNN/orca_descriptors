Интерфейс командной строки
======================

Библиотека предоставляет интерфейс командной строки (CLI) для запуска бенчмарков и оценки времени расчетов.

Команды
--------

run_benchmark
~~~~~~~~~~~~~

Запустить бенчмарк-расчет для калибровки оценки времени::

   orca_descriptors run_benchmark [ОПЦИИ]

Опции:

* ``--working_dir``: Рабочая директория для расчетов (по умолчанию: текущая директория)
* ``--functional``: DFT функционал или полуэмпирический метод (по умолчанию: AM1)
* ``--basis_set``: Базисный набор (по умолчанию: def2-SVP)
* ``--n_processors``: Количество процессоров (по умолчанию: 1)
* Все остальные параметры ORCA также доступны

Пример::

   orca_descriptors run_benchmark --working_dir ./calculations --n_processors 4

approximate_time
~~~~~~~~~~~~~~~~

Оценить время расчета для молекулы без запуска расчета::

   orca_descriptors approximate_time --molecule SMILES [ОПЦИИ]

Обязательные аргументы:

* ``--molecule``: SMILES строка молекулы

Опции:

* ``--method_type``: Тип расчета - Opt, SP или Freq (по умолчанию: Opt)
* ``--n_opt_steps``: Ожидаемое количество шагов оптимизации (для метода Opt)
* Все параметры ORCA доступны

Пример::

   orca_descriptors approximate_time \\
       --molecule CCO \\
       --method_type Opt \\
       --functional PBE0 \\
       --basis_set def2-TZVP \\
       --n_processors 8

clear
~~~~~

Удалить все файлы ORCA в рабочей директории (полезно для очистки файлов, которые не были удалены из-за ошибок)::

   orca_descriptors clear [ОПЦИИ]

Эта команда удаляет все файлы, связанные с ORCA (`.inp`, `.out`, `.log`, `.gbw`, `.cube`, `.prop` и т.д.) и все файлы, начинающиеся с ``orca_`` (включая временные файлы типа ``orca_<hash>.tmp0``, ``orca_<hash>.tmp1`` и т.д.) из рабочей директории.

Опции:

* ``--working_dir``: Рабочая директория для очистки (по умолчанию: текущая директория)
* Все остальные параметры ORCA доступны, но не используются

Пример::

   orca_descriptors clear --working_dir ./calculations

purge_cache
~~~~~~~~~~~

Удалить кеш ORCA::

   orca_descriptors purge_cache [ОПЦИИ]

Эта команда очищает директорию кеша ORCA, удаляя все закешированные результаты расчетов.

Опции:

* ``--cache_dir``: Директория кеша для очистки (по умолчанию: output_dir/.orca_cache)
* ``--output_dir``: Директория вывода (используется для определения местоположения кеша, если cache_dir не указан)
* Все остальные параметры ORCA доступны, но не используются

Пример::

   orca_descriptors purge_cache --output_dir ./calculations

Параметры удаленного кеша
--------------------------

CLI поддерживает удаленное кеширование через API. Доступен публичный сервер кеша по адресу ``https://api.orca-descriptors.massonnn.ru``.

Для использования удаленного кеша:

1. Зарегистрируйтесь на `orca-descriptors.massonnn.ru <https://orca-descriptors.massonnn.ru>`_ и выпустите API токен
2. Укажите API токен с помощью параметра ``--cache_api_token``::

   orca_descriptors run_benchmark \\
       --cache_api_token ваш_api_токен \\
       --working_dir ./calculations

Доступные параметры удаленного кеша:

* ``--cache_server_url``: URL удаленного сервера кеша (по умолчанию: https://api.orca-descriptors.massonnn.ru)
* ``--cache_api_token``: API токен для аутентификации удаленного кеша (требуется для удаленного кеша)
* ``--cache_timeout``: Таймаут для запросов удаленного кеша в секундах (по умолчанию: 30)
* ``--cache_only``: Только использовать кеш и не запускать расчеты ORCA (по умолчанию: False)

Пример с удаленным кешем::

   orca_descriptors approximate_time \\
       --molecule CCO \\
       --cache_api_token ваш_api_токен \\
       --cache_only

Доступные параметры
-------------------

Все параметры класса ``Orca`` доступны как аргументы командной строки:

* ``--script_path``: Путь к исполняемому файлу ORCA (по умолчанию: 'orca')
* ``--working_dir``: Рабочая директория для расчетов
* ``--output_dir``: Директория для выходных файлов
* ``--functional``: DFT функционал или полуэмпирический метод (по умолчанию: AM1)
* ``--basis_set``: Базисный набор (по умолчанию: def2-SVP)
* ``--method_type``: Тип расчета: Opt, SP или Freq (по умолчанию: Opt)
* ``--dispersion_correction``: Коррекция дисперсии, например, D3BJ (по умолчанию: D3BJ). Используйте 'None' для отключения.
* ``--solvation_model``: Модель сольватации, например, 'COSMO(Water)' (по умолчанию: None)
* ``--n_processors``: Количество процессоров (по умолчанию: 1)
* ``--max_scf_cycles``: Максимальное количество SCF циклов (по умолчанию: 100)
* ``--scf_convergence``: Порог сходимости SCF (по умолчанию: 1e-6)
* ``--charge``: Заряд молекулы (по умолчанию: 0)
* ``--multiplicity``: Спиновая мультиплетность (по умолчанию: 1)
* ``--cache_dir``: Директория для кеширования результатов
* ``--log_level``: Уровень логирования: DEBUG, INFO, WARNING, ERROR (по умолчанию: INFO)
* ``--max_wait``: Максимальное время ожидания создания выходного файла в секундах (по умолчанию: 300)

Справка
----

Получить справку по любой команде::

   orca_descriptors --help
   orca_descriptors run_benchmark --help
   orca_descriptors approximate_time --help
   orca_descriptors clear --help
   orca_descriptors purge_cache --help

