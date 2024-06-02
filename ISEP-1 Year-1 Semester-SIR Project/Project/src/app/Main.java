import java.io.*;
import java.util.Scanner;

import static java.lang.Thread.sleep;

public class Main {
    public static Scanner input = new Scanner(System.in);
    public static String INVALID = "Resposta inválida.%n";
    public static float[] PARAMETERS = new float[4]; //método, passo, tamanho da pop, num de dias
    public static String inputFile;
    public static String[] outputFiles; //input e outputs

    public static String[] name;
    public static float[][] RATES;
    //RATES[i][0] = beta , RATES[i][1] gama =  RATES[i][2] = ro, RATES[i][3] = alfa,

    public static void main(String[] args) throws IOException {
        mainMenu(args);
    }

    private static void mainMenu(String[] args) throws IOException {
        chooseMainMenu(args);
    }

    private static void chooseMainMenu(String[] args) throws IOException {
        if (args.length == 0) {
            interactive();
        } else if (args.length == 9) {
            boolean proceed = notInteractive(args);
            //analiseFileExtension();
            if(((PARAMETERS[0] == 1) || (PARAMETERS[0] == 2)) && (PARAMETERS[1] > 0 && PARAMETERS[1] <= 1) && PARAMETERS[2] >= 1 && PARAMETERS[3] >= 1 && proceed) {
                try {
                    int people = countNumberOfPeople();
                    if (people != 0) {
                        RATES = new float[people][4];
                        name = new String[people];
                        outputFiles = new String[people];
                        readFile();
                        for (int i = 0; i < name.length; i++) {
                            outputFiles[i] = outputFile(PARAMETERS, i);
                            program(outputFiles[i], i);
                        }
                        prepareToCreateGraphs();
                    } else {
                        System.out.println("O ficheiro inserido não tem dados suficientes para ser analisado, tente novamente");
                    }
                } catch (FileNotFoundException e) {
                    System.out.println("O ficheiro inserido não existe, verifique o seu nome.");
                }
            }else{
                System.out.println("Os valores inseridos não são válidos, tente novamente");
            }
        } else {
            System.out.println("Número de argumentos inválido, tente novamente.");
        }
    }

    public static void program(String output, int person) throws IOException {
        FileWriter outputWriter = new FileWriter(output);
        float[] y = {PARAMETERS[2] - 1, 1, 0}; // y[0] = S y[1] = I ,y[2] = R
        float[][] yn;
        switch (Math.round(PARAMETERS[0])) {
            case 1 -> {
                yn = Euler(RATES[person], y, 0, Math.round(PARAMETERS[3] / PARAMETERS[1]), PARAMETERS[1]);
                printMatrix(yn, outputWriter);
            }
            case 2 -> {
                yn = RungeKutta_4(0, y, PARAMETERS[1], Math.round(PARAMETERS[3] / PARAMETERS[1]), person);
                printMatrix(yn, outputWriter);
            }
        }
        outputWriter.close();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////

    private static float[] funcDiferencials(float[] y, float x, float[] constants) {
        float[] eachDayValue = new float[3];
        eachDayValue[0] = -constants[0] * y[0] * y[1];
        //-beta * S * I
        eachDayValue[1] = constants[2] * constants[0] * y[0] * y[1] - constants[1] * y[1] + constants[3] * y[2];
        //ro * beta * S * I - gama * I + alfa * R
        eachDayValue[2] = constants[1] * y[1] - constants[3] * y[2] + (1 - constants[2]) * constants[0] * y[0] * y[1];
        //gama * I - alfa * R + (1 - ro) * beta * S * I1
        return eachDayValue;
    }

    /**
     * Constants definition<br>
     * constants[0] = Beta, constants[1] = Alpha, constants[2] = Gama, constants[3] = Ro
     * <br><br>Variables definition<br>
     * y[0] = S, y[1] = I, s[2] = R
     *
     * @param constants the program constants
     * @param y         the y that has our 3 y
     * @param t         the x that simulates the time
     * @param n         the number of steps
     * @param h         the step
     * @return the matrix made
     */
    public static float[][] Euler(float[] constants, float[] y, float t, int n, float h) {
        float[][] resultsMatrix = new float[n][3];
        int i = 1;
        resultsMatrix[0][0] = y[0];
        resultsMatrix[0][1] = y[1];
        resultsMatrix[0][2] = y[2];
        float[] y0 = y;
        do {
            float[] temp = funcDiferencials(y0, t + i * h, constants);
            y[0] = y0[0] + h * temp[0];
            y[1] = y0[1] + h * temp[1];
            y[2] = y0[2] + h * temp[2];
            y0 = y;
            t += h;
            resultsMatrix[i][0] = y0[0];
            resultsMatrix[i][1] = y0[1];
            resultsMatrix[i][2] = y0[2];
            i++;
        } while (i < n);
        return resultsMatrix;
    }

    public static float[][] RungeKutta_4(float x0, float[] y0, float h, int n, int person) {
        y0[0] = PARAMETERS[2] - 1;
        y0[1] = 1;
        y0[2] = 0;
        float[][] resultsMatrix = new float[n][3];
        int j = 1;
        resultsMatrix[0][0] = y0[0];
        resultsMatrix[0][1] = y0[1];
        resultsMatrix[0][2] = y0[2];
        float[] y = new float[3];

        do {
            float[] k1 = calculateK1(h, y0, x0, person);
            float[] k2 = calculateK2(h, y0, x0, k1, person);
            float[] k3 = calculateK3(h, y0, x0, k2, person);
            float[] k4 = calculateK4(h, y0, x0, k3, person);
            float[] k = calculateK(k1, k2, k3, k4);
            y[0] = y0[0] + k[0];
            y[1] = y0[1] + k[1];
            y[2] = y0[2] + k[2];
            y0 = y;
            resultsMatrix[j][0] = y0[0];
            resultsMatrix[j][1] = y0[1];
            resultsMatrix[j][2] = y0[2];
            x0 += h;
            j++;
        } while (j < n);
        return resultsMatrix;
    }

    private static float[] calculateK1(float h, float[] y0, float x0, int person) {
        float[] k1 = new float[3];
        float[] temp = funcDiferencials(y0, x0, RATES[person]);
        k1[0] = h * temp[0];
        k1[1] = h * temp[1];
        k1[2] = h * temp[2];
        return k1;
    }

    private static float[] calculateK2(float h, float[] y0, float x0, float[] k1, int person) {
        float[] k2 = new float[3];
        float[] y = new float[3];
        y[0] = y0[0] + k1[0]/2;
        y[1] = y0[1] + k1[1]/2;
        y[2] = y0[2] + k1[2]/2;
        float[] temp = funcDiferencials(y, x0+h/2, RATES[person]);
        k2[0] = h * temp[0];
        k2[1] = h * temp[1];
        k2[2] = h * temp[2];
        return k2;
    }

    private static float[] calculateK3(float h, float[] y0, float x0, float[] k2, int person) {
        float[] k3 = new float[3];
        float[] y = new float[3];
        y[0] = y0[0] + k2[0]/2;
        y[1] = y0[1] + k2[1]/2;
        y[2] = y0[2] + k2[2]/2;
        float[] temp = funcDiferencials(y, x0+h/2, RATES[person]);
        k2[0] = h * temp[0];
        k2[1] = h * temp[1];
        k2[2] = h * temp[2];
        return k3;
    }

    private static float[] calculateK4(float h, float[] y0, float x0, float[] k3, int person) {
        float[] k4 = new float[3];
        float[] y = new float[3];
        y[0] = y0[0] + k3[0];
        y[1] = y0[1] + k3[1];
        y[2] = y0[2] + k3[2];
        float[] temp = funcDiferencials(y, x0+h/2, RATES[person]);
        k3[0] = h * temp[0];
        k3[1] = h * temp[1];
        k3[2] = h * temp[2];
        return k4;
    }

    private static float[] calculateK(float[] k1, float[] k2, float[] k3, float[] k4) {
        float[] k = new float[3];
        k[0] = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6;
        k[1] = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6;
        k[2] = (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6;
        return k;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////

    private static boolean notInteractive(String[] args) {
        inputFile = args[0];
        boolean valid = true;
        for (int i = 2; i < args.length && valid; i += 2) {
            switch (args[i - 1]) {
                case "-m" -> PARAMETERS[0] = Float.parseFloat(args[i]); //método: 1 - Euler; 2 - Runge Kutta
                case "-p" -> PARAMETERS[1] = Float.parseFloat(args[i]); //passo
                case "-t" -> PARAMETERS[2] = Float.parseFloat(args[i]); //tamanho da população
                case "-d" -> PARAMETERS[3] = Float.parseFloat(args[i]); //número de dias
                default -> {
                    System.out.println("Argumento inserido foi inválido");
                    valid = false;
                }
            }
        }
        return valid;
    }

    public static void interactive() throws IOException {
        interactiveMenu();
    }

    private static void interactiveMenu() throws IOException { // menu de início (interativo)
        boolean exit = false;
        while (!exit) {
            displayInteractiveMenu();
            exit = chooseInteractiveMenu(); // mudar de menu
        }
    }

    private static void displayInteractiveMenu() { //incompleto
        System.out.printf("Cálculo da propagação de fake news perante as condições que a influenciam.%n%n");
        System.out.printf("Introduza o número correspondente ao menu a aceder.%n");
        System.out.printf("(1) - Introduzir parâmetros.%n");
        System.out.printf("(2) - Ver parâmetros.%n");
        System.out.printf("(3) - Executar programa%n");
        System.out.printf("(4) - Sair.%n");
    }

    public static void printMatrix(float[][] resultsMatrix, FileWriter outputWriter) throws IOException {
        StringBuilder phrase = new StringBuilder("");
        phrase.append("Dia;S;I;R;N\n");
        float total;

        for (int i = 0; i < resultsMatrix.length; i++) {
            if (i % (1/PARAMETERS[1]) == 0) {
                total = resultsMatrix[i][0] + resultsMatrix[i][1] + resultsMatrix[i][2];
                phrase.append((int) (i * PARAMETERS[1])).append(";").append(resultsMatrix[i][0]).append(";").append(resultsMatrix[i][1]).append(";").append(resultsMatrix[i][2]).append(";").append(total).append("\n");
            }
        }
        outputWriter.write(String.valueOf(phrase));
    }

    private static boolean chooseInteractiveMenu() throws IOException {
        int menu = scanInt(); // escolher menu
        boolean validAnswer = false;
        cls();
        while (!validAnswer) {
            switch (menu) {
                case 1: // entrar no menu de introdução de parâmetros
                    printInputParameteres();
                    chooseParameteresMenu();
                    validAnswer = true;
                    break;
                case 2: // print ao Array de parametros
                    printParameters();
                    validAnswer = true;
                    break;
                case 3:
                    if(PARAMETERS[0] != 0 && PARAMETERS[1] != 0 && PARAMETERS[2] != 0 && PARAMETERS[3] != 0 && inputFile != null && !inputFile.isEmpty()) {
                        for (int i = 0; i < name.length; i++) {
                            outputFiles[i] = outputFile(PARAMETERS, i);
                            program(outputFiles[i], i);
                        }
                        prepareToCreateGraphs();
                    }else{
                        System.out.println("Ainda não inseriu valores suficientes para executar o programa, insira primeiro os valores necessários");
                        exit();
                        validAnswer = true;
                    }
                    break;
                case 4: // Sair do programa
                    return true;
                default:
                    System.out.printf(INVALID);
                    exit();
                    displayInteractiveMenu();
                    menu = scanInt();
                    break;
            }
        }
        return false;
    }

    public static String outputFile(float[] parameters, int person) {
        String outputname = name[person];
        outputname += "m";
        outputname += String.valueOf((int) parameters[0]);
        outputname += "p";
        outputname += String.valueOf(parameters[1]);
        outputname += "t";
        outputname += String.valueOf((int) parameters[2]);
        outputname += "d";
        outputname += String.valueOf((int) parameters[3]);
        outputname = outputname.replace(".", "");
        outputname += ".csv";
        return (outputname);
    }

    private static int countNumberOfPeople() throws FileNotFoundException {
        File csvFile = new File(inputFile);
        Scanner sc = new Scanner(csvFile);
        sc.nextLine();
        int i = 0;
        while(sc.hasNext()){
            sc.nextLine();
            i++;
        }
        sc.close();
        return i;
    }

    private static void readFile() throws FileNotFoundException {
        File csvFile = new File(inputFile);
        Scanner sc = new Scanner(csvFile);
        sc.nextLine();
        for (int i = 0; sc.hasNext(); i++) {
            String values = sc.nextLine();
            String[] valuesFromLine = values.split(";");
            name[i] = valuesFromLine[0];
            for (int j = 0; j < RATES[i].length; j++) {
                valuesFromLine[j + 1] = valuesFromLine[j + 1].replace(",", ".");
                RATES[i][j] = Float.parseFloat(valuesFromLine[j + 1]);
            }
        }
        sc.close();
    }

    private static void printInputParameteres() {
        System.out.printf("Selecionou o menu: (1) - Introduzir parâmetros.%n");
        System.out.printf("Introduza o número correspondente ao parâmetro a introduzir.%n");
        System.out.printf("(1) - Escolher método.%n");
        System.out.printf("(2) - Definir passo.%n");
        System.out.printf("(3) - Definir tamanho da população.%n");
        System.out.printf("(4) - Definir número de dias.%n");
        System.out.printf("(5) - Escolher o ficheiro de input com taxas.%n");
        System.out.printf("(6) - Voltar ao início.%n");
    }

    private static void chooseParameteresMenu(){
        int parameter = scanInt(); // escolher parâmetro
        cls();
        while (parameter != 6) {
            switch (parameter) {
                case 1 -> {
                    System.out.printf("Selecionou o parâmetro: (1) - Método.%n");
                    PARAMETERS[parameter - 1] = (float) scanMethod();
                    printInputParameteres();
                    parameter = scanInt();
                    cls();
                }
                case 2 -> {
                    System.out.printf("Selecionou o menu: (2) - Definir passo.%n");
                    PARAMETERS[parameter - 1] = scanFloat();
                    printInputParameteres();
                    parameter = scanInt();
                    cls();
                }
                case 3 -> {
                    System.out.printf("Selecionou o parâmetro: (3) - Tamanho da População%n");
                    printValue();
                    PARAMETERS[parameter - 1] = (float) scanInt();
                    printInputParameteres();
                    parameter = scanInt();
                    cls();
                }
                case 4 -> {
                    System.out.printf("Selecionou o parâmetro: (4) - Número de dias%n");
                    printValue();
                    PARAMETERS[parameter - 1] = (float) scanInt();
                    printInputParameteres();
                    parameter = scanInt();
                    cls();
                }
                case 5 -> {
                    System.out.printf("Selecionou o parâmetro: (5) - Ficheiro de input%n");
                    inputFile = scanString();
                    //analiseFileExtension();
                    try {
                        int people = countNumberOfPeople();
                        if (people != 0) {
                            RATES = new float[people][4];
                            name = new String[people];
                            outputFiles = new String[people];
                            readFile();
                            printInputParameteres();
                            parameter = scanInt();
                            cls();
                        } else {
                            System.out.println("O ficheiro inserido não tem dados suficientes para ser analisado, tente novamente");
                        }
                    } catch (FileNotFoundException e) {
                        System.out.println("O ficheiro inserido não existe, tente novamente");
                    }
                }
                default -> {
                    System.out.printf(INVALID);
                    exit();
                    printInputParameteres();
                    parameter = scanInt();
                }
            }
        }
    }

    private static int scanMethod() {
        int n = 0;
        boolean valido;
        System.out.printf("Introduza o número \"1\" para o método de Euler.%n");
        System.out.printf("Introduza o número \"2\" para o método de Runge Kutta de 4ª ordem.%n");
        System.out.printf("Não poderá avançar até introduzir um valor válido.%n");
        do {
            valido = true;
            try {
                n = Integer.parseInt(input.nextLine());
            } catch (NumberFormatException e) {
                valido = false;
                System.out.printf(INVALID);
            }
            if (valido && (n != 1) && (n != 2)){
                valido = false;
                System.out.printf(INVALID);
            }
        } while (!valido);
        cls();
        return n;
    }

    private static void printParameters() {
        String metodo = "Nulo.";
        if (PARAMETERS[0] == 1) {
            metodo = "Euler";
        }
        if (PARAMETERS[0] == 2) {
            metodo = "Runge Kutta de 4ª ordem";
        }
        System.out.printf("Selecionou o menu: (2) - Ver parâmetros.%n");

        System.out.printf("Método: %s%nPasso: %f%nPopulação: %f%nDias: %f%nNome do ficheiro de input: %s%n%n", metodo, PARAMETERS[1], PARAMETERS[2], PARAMETERS[3], inputFile);
        if(inputFile != null && !inputFile.isEmpty()){
            String b = "Beta", g = "Gama", r = "Ró", a = "Alfa", empty = "";

            System.out.printf("%-8s|   %-6s |   %-6s |    %-5s |   %-6s%n", empty, b, g, r, a);
            System.out.println("----------------------------------------------------");
            for(int i = 0; i < name.length ; i++) {
                System.out.printf("%-7s | %-8f | %-8f | %-8f | %-8f%n", name[i], RATES[i][0], RATES[i][1], RATES[i][2], RATES[i][3]);
            }
        }else{
            System.out.println("Os valores de taxas ainda não podem ser visualizadas");
        }
        System.out.println();
        exit();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    private static void exit() {
        System.out.printf("Use a tecla ENTER para voltar atrás.%n");
        input.nextLine();
        cls();
    }

    private static void printValue() {
        System.out.printf("Introduza um valor.%n");
    }

    private static void cls() {
        for (int i = 0; i < 50; i++) {
            System.out.printf("%n");
        }
    }

    private static int scanInt() { // dar Scan a um inteiro pelo teclado
        int intNumber = 0;
        boolean valido;
        do {
            valido = true;
            try {
                intNumber = Integer.parseInt(input.nextLine());
                if(intNumber <= 0){
                    valido = false;
                    System.out.println("Valor inválido. Insira um valor válido.");
                }
            } catch (NumberFormatException e) {
                valido = false;
                System.out.println("Valor inválido. Insira um valor válido.");
            }
        }while(!valido);
        cls();
        return intNumber;
    }

    private static float scanFloat() { // dar Scan a um double pelo teclado
        float floatNumber = 0;
        boolean valido;
        printValue();
        System.out.printf("Não poderá avançar até introduzir um valor válido, pertencente ao intervalo ]0,1].%n");
        do {
            valido = true;
            try {
                floatNumber = Float.parseFloat(input.nextLine());
            } catch (NumberFormatException e) {
                valido = false;
                System.out.printf(INVALID);
            }
            if (valido && (floatNumber <= 0 || floatNumber > 1)){
                valido = false;
                System.out.printf(INVALID);
            }
        } while (!valido);
        cls();
        return floatNumber;
    }

    private static String scanString() {
        String Word;
        Word = input.next();
        input.nextLine();
        cls();
        return Word;
    }

    private static void prepareToCreateGraphs() throws IOException {
        try {
            sleep(10000);
        } catch (InterruptedException e) {
            System.out.println("no sleep");
        }
        for (int i = 0; i < name.length; i++) {
            createGraph(outputFiles[i], "Graph" + name[i] + ".png");
        }
    }
    private static void createGraph(String file, String output) throws IOException{
        ProcessBuilder procedure = new ProcessBuilder("gnuplot", "-p");
        Process activity = procedure.start();

        BufferedWriter gnuplot = new BufferedWriter(new OutputStreamWriter(activity.getOutputStream()));
        BufferedReader errorReader = new BufferedReader(new InputStreamReader(activity.getErrorStream()));

        gnuplot.write("set xtics (`stats '" + file + "' using 1 skip 1 nooutput`)\n");
        gnuplot.write("set ytics (`stats '"+ file +"' using 5 skip 1 nooutput`)\n");
        gnuplot.write("set datafile separator ';'\n");
        gnuplot.write("set xlabel \"Dias\"\n");
        gnuplot.write("set ylabel \"Pessoas\"\n");
        gnuplot.write("set autoscale ';'\n");
        gnuplot.write("set grid ';'\n");
        gnuplot.write("set terminal png\n");
        gnuplot.write("seset output '" + output + "'\n");
        gnuplot.write("plot '"+ file +"' using 1:2 skip 1 lw 3 lc 7 with lines title 'S', '' using 1:3 skip 1 lw 3 lc 5 with lines title 'I', '' using 1:4 skip 1 lw 3 lc 2 with lines title 'R'\n");

        errorReader.close();
        gnuplot.close();
    }
}