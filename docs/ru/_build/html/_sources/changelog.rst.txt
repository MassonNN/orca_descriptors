История изменений
=========


Версия 0.2.2
--------------

Добавлено
~~~~~

* Добавлено ``numpy>=1.20.0`` to project dependencies (numpy was used but not declared)
* Добавлено dynamic time estimation updates in batch processing - time estimates are now refined based on actual execution times of previous molecules
* Добавлено ``_get_available_descriptors()`` method to dynamically discover available descriptor methods

Изменено
~~~~~~~

* Обновлен парсер дипольного момента для приоритизации газофазных значений при их наличии (для расчетов без сольватации)
* Улучшен алгоритм оценки времени:
  * Изменено scaling exponent from O(N^3.5) to O(N^2.5) for more realistic estimates
  * Используется ``total_time`` из benchmark вместо ``scf_time`` как базовая единица
  * Более реалистичная оценка шагов оптимизации (15-35 шагов вместо 10-50)
  * Удален искусственный предел в 24 часа
* Рефакторинг метода ``calculate_descriptors()``:
  * Убрано дублирование кода (заменена большая цепочка if-elif на вызовы методов через ``getattr``)
  * Удален избыточный список ``all_descriptors`` - дескрипторы теперь обнаруживаются динамически
  * Удалены ненужные комментарии
  * Улучшена поддерживаемость и читаемость кода

Исправлено
~~~~~

* Исправлено dipole moment parser to correctly extract gas-phase values from ORCA output when available
* Исправлено time estimation showing unrealistic values (e.g., 47 hours for 2 molecules) - now provides accurate estimates based on actual benchmark data

Технические детали
~~~~~~~~~~~~~~~~~

* Оценщик времени теперь использует экспоненциальное скользящее среднее для лучшей точности прогнозирования
* Методы дескрипторов вызываются динамически с использованием ``getattr(self, desc_name)``
* Автоматическое обнаружение дескрипторов устраняет необходимость поддерживать списки дескрипторов вручную

