Установка
=========

Требования
----------

* Python >= 3.10
* ORCA 6.0.1 установлен и доступен в PATH
* RDKit >= 2023.0.0

Установка ORCA
--------------

ORCA — это программа квантовой химии, которую необходимо установить отдельно.

1. Скачайте ORCA с официального сайта: https://orcaforum.kofo.mpg.de/
2. Распакуйте и добавьте исполняемый файл ORCA в системный PATH
3. Проверьте установку, выполнив::

   orca --version

Установка библиотеки
---------------------

С помощью pip
~~~~~~~~~~~~~

.. code-block:: bash

   pip install orca-descriptors

С помощью Poetry (для разработки)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone <repository-url>
   cd orca_descriptors
   poetry install

Проверка
--------

После установки проверьте, что библиотека работает::

   python -c "from orca_descriptors import Orca; print('OK')"

Также можно проверить CLI::

   orca_descriptors --help

