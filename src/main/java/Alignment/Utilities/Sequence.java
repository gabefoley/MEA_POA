package Alignment.Utilities;

/**
 * Created by gabe on 7/09/2016.
 */
public class Sequence {

    private String name;
    private String content;

    public Sequence(String name, String content){
        this.name = name;
        this.content = content;
    }

    public String getName(){
        return this.name;
    }


    public String getContent(){
        return  this.content;
    }

}
