# Projekt: Rozwiązywanie układów liniowych

## Opis

Projekt ten zawiera implementacje trzech różnych metod rozwiązywania układów liniowych: metody Jacobiego, metody Gaussa-Seidla oraz faktoryzacji LU. Układy liniowe są rozwiązywane dla macierzy kwadratowej \( A \) o rozmiarze \( N \) oraz wektora \( b \). 

## Metody

### 1. Metoda Jacobiego

Metoda Jacobiego jest iteracyjną metodą rozwiązywania układów równań liniowych. Polega na iteracyjnym poprawianiu przybliżenia rozwiązania na podstawie poprzedniego przybliżenia.

### 2. Metoda Gaussa-Seidla

Metoda Gaussa-Seidla również jest iteracyjną metodą rozwiązywania układów równań liniowych. Podobnie jak metoda Jacobiego, polega na iteracyjnym poprawianiu przybliżenia rozwiązania, jednakże wykorzystuje najnowsze dostępne wartości w obliczeniach.

### 3. Faktoryzacja LU

Faktoryzacja LU polega na podziale macierzy \( A \) na iloczyn dwóch macierzy \( L \) i \( U \), gdzie \( L \) to macierz trójkątna dolna, a \( U \) to macierz trójkątna górna. Rozwiązanie układu równań sprowadza się do rozwiązania dwóch układów równań trójkątnych.

## Plik projektowy

- **projekt2.cpp**: Plik zawiera implementacje powyższych metod oraz kod główny programu, który testuje działanie metod i zbiera wyniki.

## Uruchamianie

Aby uruchomić program, należy skompilować plik `projekt2.cpp` z wykorzystaniem kompilatora obsługującego standard C++. Następnie można uruchomić skompilowany plik wykonywalny.

Przykładowe polecenie kompilacji za pomocą g++:

```bash
g++ projekt2.cpp -o projekt2
```

Aby uruchomić skompilowany program:

```bash
./projekt2
```

## Wyniki

- Wyniki działania programu zostaną wyświetlone na standardowym wyjściu oraz zapisane do plików tekstowych.
- Dla metod iteracyjnych (Jacobi i Gauss-Seidel) obliczany jest czas wykonania oraz norma residuum dla różnych wartości parametrów.
- Dla metody faktoryzacji LU obliczany jest czas wykonania oraz norma residuum.

## Wymagania

- Kompilator obsługujący standard C++.
- Biblioteki standardowe C++: iostream, cmath, chrono, fstream.

## Autor

Autor: [Paulina Machcińska](https://github.com/PaulinaM122)
