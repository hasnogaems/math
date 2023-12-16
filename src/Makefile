CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Werror
CCOVTEST=-fprofile-arcs -ftest-coverage
LIBS = -lm

all: s21_math.a test gcov_report

s21_math.a: s21_math.o
	ar rcs $@ $^

s21_math.o: s21_math.c s21_math.h
	$(CC) $(CFLAGS) -c $< -o $@

test: test_s21_math.o s21_math.o
	$(CC) $(CFLAGS) $(CCOVTEST) $^ $(LIBS) -o $@ -lcheck -lm -lpthread
	./test

test_s21_math.o: test_s21_math.c s21_math.h
	$(CC) $(CFLAGS) $(CCOVTEST) -I. -c $< -o $@

gcov_report: test s21_math.a
	$(CC) $(CFLAGS) -pedantic -lcheck $(CCOVTEST) test_s21_math.c s21_math.a -o $@
	./gcov_report
	lcov -t gcov_report -o test.info -c -d .
	genhtml -o report test.info
	open report/index.html

clang:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -n *.c *.h
	
clean:
	rm -f *.a
	rm -f *.o
	rm -f test
	rm -f *.gc*
	rm -f gcov_report.html
	rm -f gcov_report
	rm -f test.info
	rm -f -R report
rebuild: clean all
	clear
