1. Zmiany w README:
- wersje wymaganych narzędzi zewnętrznych np. GATK i Annovara  ✔
- przykładowe wywołanie
2. Dodatkowe funkcjonalności:
- policzyć korelację między genotypem a sygnałem H3K27ac w obrębie promotora lub enhancera, w którym leży SNP ✔
- umożliwić uzyskanie częściowych wyników jeśli użytkownik ma tylko dane o wariantach ✔
3. Zmiany w formacie tabeli wynikowej:
- kolumny:
CHROM, POS, REF, ALT, AC, AF, AN, Num homref, Num het, Num homalt, gnomAD_genome_ALL, gnomad_for_selected_population, binom_pval, corrected_binom_pval, genomic element, Gene_name, Gene_ID, Transcript_ID, p-value_for_correlation_genotype_vs_expression, p-value_for_correlation_genotype_vs_h3k27ac, p-value_for_correlation_expression_vs_h3k27ac(enhancers only)  ✔
- p-wartości pokazane do 3 miejsc po przecinku ✔
- rozdzielić wyniki dla różnych transkryptów na osobne wiersze ✔
4. Zmiany na wykresach:
- zapis wariantu: chr1:123456A>T ✔
- etykiety na osi X najlepiej jakby odpowiadały genotypom np. AA AT TT ✔
- dodać wykresy dla korelacji genotypu z H3K27ac ✔
5. Zapisać w pliku domyślne przypisania genów do enhancerów, żeby nie przeliczać ich ciągle na nowo ✔


Pytania:
- czy połączyć tabelę outputową enhancerów i promotorów --można
