package CrazyE;

import CrazyE.graphic.CEMinimapFragment;
import CrazyE.world.environment.*;
import CrazyE.graphic.CEMinimapRenderer;
import mindustry.core.*;
import mindustry.mod.*;

import java.lang.reflect.*;

import static mindustry.io.SaveVersion.*;
import static mindustry.Vars.*;

public class CrazyEnvironment extends Mod{
    public static CustomChunk Gravity;

    public CrazyEnvironment(){}

    @Override
    public void loadContent(){
        Gravity = new Gravity();
    }

    @Override
    public void init() {
        super.init();

        try{
            Field f = Renderer.class.getDeclaredField("minimap");
            f.setAccessible(true);
            f.set(renderer, new CEMinimapRenderer());
            ui.minimapfrag = new CEMinimapFragment();
            ui.minimapfrag.build(ui.hudGroup);
        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }
}

