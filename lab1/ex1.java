import java.util.Scanner;

public class ex1 {
    public static String ex1(String sequence) {
        StringBuilder alphabet = new StringBuilder();
        for (char c : sequence.toCharArray()) {
            if (Character.isLetter(c) && alphabet.indexOf(String.valueOf(c)) == -1) {
                alphabet.append(c);
            }
        }
        return alphabet.toString();
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        System.out.print("sequnce: ");
        String input = scanner.nextLine();
        String result = ex1(input);
        System.out.println("alphabet: " + result);
        scanner.close();

    }
}