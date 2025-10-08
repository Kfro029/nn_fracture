#By Alena V. Favorskaya. 25.07.2025

сборка rect.
при копировании от человека исходников, папка build удаляется. далее
cd rect
ln -s ../advlib
ln -s ../mlib/mlib
./build.sh

заменить в файле run.sh
/home/aleanera/work/rect/build/rect
на актуальный путь до исполняемого файла
/rect/build/rect

затем запуск
chmod +x run.sh
./run.sh



**** для себя
--Шаг 1. Собрать rect
cd rect
ln -s ../advlib
ln -s ../mlib/mlib
./build.sh


--Шаг 2. Сгенерировать трещены
--generate_random.sh генерирует $TOTAL_FILES файлов, в которых рандомно от 1 до 10 трещин. --Для его запуска:
chmod +x generate_random.sh
./generate_random.sh

--Шаг 3. Запуск расчета
-- chmod +x run.sh
-- ./run.sh


--для запуска на сервере
--подключение на сервер
ssh login@ip_адрес_сервера
--перенос файлов на серевер
scp fracture_generator_vtk.py login@ip_адрес_сервера:~/путь/на/сервере/
scp conf.conf login@ip_адрес_сервера:~/путь/на/сервере/
--запуск

--скачивание на ПК
scp login@ip_адрес_сервера:/путь/к/fractures.vtk ~/Загрузки/
