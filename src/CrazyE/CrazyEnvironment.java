package CrazyE;

import CrazyE.world.DynamicEnvironment.*;
import mindustry.mod.*;
import static mindustry.io.SaveVersion.*;

public class CrazyEnvironment extends Mod{
    public static CustomChunk Gravity;

    public CrazyEnvironment(){

    }

    @Override
    public void loadContent(){
        Gravity = new Gravity();
    }
}
