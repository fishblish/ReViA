# search-for-regulatory-variants
Przygotowywanie narzędzia do poszukiwania wariantów regulatorowych w genomie.

Wszystkie kluczowe funkcje znajdują się w pliku firstmodule.py, wywołanie ich odbywa się w pliku main.py.

Kolejne kroki występujące w analizie:
- Wybór wariantów, które znajdują się w regionach regulatorowych i są bialleliczne.
- Przypisanie częstości występowania wariantu w populacji i przefiltrowanie po tej wartości. Jest to wykonane za pomoca ANNOVARu.
- Wykonanie testu dwumianowego - wybór wariantów, których jest istotnie więcej w grupie badanej niż w populacji. Wykonanie korekty Benjaminiego-Hochberga.
- Przypisanie genów do promotorów i enhancerów intronowych.
- Przypisanie genów do enhancerów: najbliższych, będących w kontakcie chromatynowym.
- Wyznaczenie korelacji między aktywnością enhancera (sygnał h3k27ac), a ekspresją genu. Wybór najlepszego genu.
- Wyznaczenie korelacji między genotypem, a ekspresją genu (na poziomie transkryptu, następnie wybór najlepszego transkryptu).
  




