# bio-alignment
Реализация алгоритмов биоинформатики Нидлмана-Вунша, Смита-Ватермана, Хиршберга и аффинного алгоритма для выравнивания биологических последовательностей.

### Технологии

Алгоритм написан на чистом `C89` без каких-либо зависимостей.

### Установка

Для `bio-alignment` требуется только `C89`-совместимый компилятор и утилита` make`.

```sh
$ cd bio-alignment
$ make
$ ./bin/bio-alignment --help
$ ./bin/bio-alignment -i data/in.fasta -o out.fasta -s blosum62 -g -5
```

### Тесты

В папке `data` расположены файлы `in1.fasta`, `in2.fasta`, `in1_1.fasta`, `in1_2.fasta`, `in2_1.fasta`, `in2_2.fasta`используемые для тестирования утилиты.

Тест `in.fasta` может быть запущен следующим образом:
```sh
$ ./bin/bio-alignment -i data/in1.fasta
```

Тест `in1.fasta`, `in2.fasta` может быть запущен следующим образом:
```sh
$ ./bin/bio-alignment -i data/in1_1.fasta,data/in1_2.fasta
```
