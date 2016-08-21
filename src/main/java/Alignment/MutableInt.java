package Alignment;

/**
 * Created by gabe on 16/08/2016.
 */
public class MutableInt {

    int value = 1;
    public void increment() {
        ++value;
    }

    public int getValue() {
        return value;
    }

    public String toString(){
        return Integer.toString(value);
    }
}
