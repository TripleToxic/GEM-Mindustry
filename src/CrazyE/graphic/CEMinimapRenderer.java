package CrazyE.graphic;

import arc.*;
import arc.graphics.*;
import arc.struct.*;
import arc.util.*;
import mindustry.content.*;
import mindustry.game.EventType.*;
import mindustry.graphics.MinimapRenderer;
import mindustry.world.*;
import CrazyE.graphic.CEMinimapField.*;

import static mindustry.Vars.*;

public class CEMinimapRenderer extends MinimapRenderer{
    private static final float updateInterval = 2f;

    private float updateCounter = 0f;
    private CEMinimapField current;
    public boolean normalMinimap, isReset;
    public FieldType field = FieldType.Heatmap;
    public String currentMap = "normal";

    public static final ObjectMap<String, CEMinimapField> fields = new ObjectMap<>();

    public CEMinimapRenderer(){
        Events.on(WorldLoadEvent.class, event -> {
            reset();
            updateAll();
            changed(currentMap);
        });

        Events.on(TileChangeEvent.class, event -> {
            if(!ui.editor.isShown() && normalMinimap){
                update(event.tile);

                //update floor below block.
                if(event.tile.block().solid && event.tile.y > 0 && event.tile.isCenter()){
                    event.tile.getLinkedTiles(t -> {
                        Tile tile = world.tile(t.x, t.y - 1);
                        if(tile != null && tile.block() == Blocks.air){
                            update(tile);
                        }
                    });
                }
            }
        });

        Events.on(TilePreChangeEvent.class, e -> {
            //update floor below a *recently removed* block.
            if(e.tile.block().solid && e.tile.y > 0 && normalMinimap){
                Tile tile = world.tile(e.tile.x, e.tile.y - 1);
                if(tile.block() == Blocks.air){
                    Time.run(0f, () -> update(tile));
                }
            }
        });

        Events.on(BuildTeamChangeEvent.class, event -> {
            if(normalMinimap) update(event.build.tile);
        });

        addField("normal", null);
    }

    public static void addField(String s, CEMinimapField f){
        fields.put(s, f);
    }

    public void changed(String s){
        currentMap = s;
        if(currentMap.equals("normal") && isReset){
            normalMinimap = true;
            isReset = false;
            return;
        }
        if(fields.containsKey(currentMap)){
            current = fields.get(currentMap);
            updateCounter = 0;
            normalMinimap = false;
        }
    }

    @Override
    public void update(){
        if((updateCounter += Time.delta) >= updateInterval){
            updateCounter %= updateInterval;

            if(!isReset){
                super.updateAll();
                isReset = true;
                return;
            }

            if(normalMinimap){
                super.update();
                super.update(); // to complete the 2 counter cycle from the super method
                return;
            }

            Pixmap p = getPixmap();

            switch(field){
                case Heatmap -> {
                    int color;
                    for(int i=0; i<world.height(); i++){
                        for(int j=0; j<world.width(); j++){
                            color = current.color(i, j);
                            p.set(j, p.height - i - 1, color);
                        }
                    }
                    getTexture().draw(p);
                }

                /*case Vectorfield -> {

                }*/
            }
        }
    }

    @Override
    public void reset(){
        super.reset();
        currentMap = "Gravity";
        normalMinimap = true;
        isReset = false;
    }
}
