package Alignment.Utilities;

/**
 * Created by gabe on 23/09/2016.
 */
public class Scratch {

    public static void main(String[] args){

        int a = 70;
        int[] b = new int[5];
        Double c = 0.0;
//        System.out.println("a is: " + a);
//        System.out.println("b is: " + b[1]);
//        System.out.println("c is: " + c);
//
//        change(a, b, c);
//
//        System.out.println("a is: " + a);
//        System.out.println("b is: " + b[1]);
//        System.out.println("c is: " + c);

        System.out.println("loadData(\n" +
                "    {\n" +
                "        name: 'graph1',\n" +
                "        nodes: [");
        makeGraph(2000);
        System.out.println("        \n" +
                "        ],\n" +
                "        links: [");
        makeLinks(2000);

        System.out.println("        ]\n" +
                "    }\n" +
                ");");

    }

    static void change(int a, int[] b, Double c){
        a = 72;
        b[1] = 63;
        c = 10.23;
    }

    static void makeGraph(int count){
        for (int i = 0; i < count; i++){
            System.out.println("{ id: 'node" + i + "', value: { label: 'node" + i + "' } },");

        }
    }

    static void makeLinks(int count){
        for (int i = 0; i < count - 1; i++){
            System.out.println("{ u: 'node" + i + "', v: 'node" + (i+1) + "', value: { label: 'link" + i + "' } },");

        }
    }
}
