package CrazyE.world.environment;

import CrazyE.graphic.*;
import arc.Events;
import arc.func.Cons;
import arc.math.Mathf;
import arc.math.geom.*;
import arc.struct.Seq;
import arc.util.*;
import mindustry.game.*;
import mindustry.gen.Unit;
import mindustry.io.SaveFileReader.*;

import java.io.*;

import static java.lang.Math.*;
import static mindustry.Vars.*;
import static mindustry.game.EventType.*;
import static mindustry.io.SaveVersion.addCustomChunk;
import static CrazyE.graphic.CEMinimapRenderer.*;

public class Gravity implements CustomChunk, CEMinimapField {
    private static final double cons = cbrt(2),
            w1 = -(float)(cons/(2d - cons)),
            w2 = (float)(1d/(2d - cons));
    private static final float c = 15f * tilesize * tilesize / 60f, maxSpeed = 0.9f*c, csquared = c * c, G = 0.1f, k = 2f, stabilizer = 0.02f, timescale = 1/64f,
    c1 = (float)(0.5d*w2),
    c2 = (float)(0.5d*(w1+w2)),
    c3 = c2,
    c4 = c1,

    d1 = (float)w2,
    d2 = (float)w1,
    d3 = (float)w2
    ;

    private static final float drawRadius = 4f;
    private int sizeX, sizeY;
    private GEMgrid PotentialField, vPotentialField, BufferField, CurrentField, EMGravityField;
    private boolean init = false;

    String debug = "";

    public Gravity(){
        Events.on(WorldLoadEvent.class, e -> {
            this.sizeX = world.width();
            this.sizeY = world.height();

            PotentialField = new GEMgrid(sizeY, sizeX);
            vPotentialField = new GEMgrid(sizeY, sizeX);
            BufferField = new GEMgrid(sizeY, sizeX);
            CurrentField = new GEMgrid(sizeY, sizeX);
            EMGravityField = new GEMgrid(sizeY, sizeX);
            init = true;
        });

        Events.run(Trigger.afterGameUpdate, this::update);

        addField("Gravity", this);
        addCustomChunk("Gravity", this);
    }

    public void EMGravity(){
        float Ex, Ey, B;
        for(int y=1; y<sizeY-1; y++){
            for(int x=1; x<sizeX-1; x++){
                Ex = -0.5f*(PotentialField.get(y, x + 1).x - PotentialField.get(y, x - 1).x);
                Ey = -0.5f*(PotentialField.get(y + 1, x).x - PotentialField.get(y - 1, x).x);

                Ex -= vPotentialField.get(y, x).y/c;
                Ey -= vPotentialField.get(y, x).z/c;

                B = 0.5f*(PotentialField.get(y, x + 1).z - PotentialField.get(y, x - 1).z)
                        - 0.5f*(PotentialField.get(y + 1, x).y - PotentialField.get(y - 1, x).y);
                EMGravityField.get(y, x).set(Ex, Ey, B);
            }
        }

        for(int x=1; x<sizeX-1; x++){
            Ex = -0.5f*(PotentialField.get(0, x + 1).x - PotentialField.get(0, x - 1).x);
            Ey = -vPotentialField.get(0, x).x;

            Ex -= vPotentialField.get(0, x).y/c;
            Ey -= vPotentialField.get(0, x).z;

            B = 0.5f*(PotentialField.get(0, x + 1).z - PotentialField.get(0, x - 1).z)
                    - vPotentialField.get(0, x).y*c;

            EMGravityField.get(0, x).set(Ex, Ey/c, B);

            Ex = -0.5f*(PotentialField.get(sizeY-1, x + 1).x - PotentialField.get(sizeY-1, x - 1).x);
            Ey = vPotentialField.get(sizeY-1, x).x;

            Ex -= vPotentialField.get(sizeY-1, x).y/c;
            Ey -= -vPotentialField.get(sizeY-1, x).z;

            B = 0.5f*(PotentialField.get(sizeY-1, x + 1).z - PotentialField.get(sizeY-1, x - 1).z)
                     + vPotentialField.get(sizeY-1, x).z/c;

            EMGravityField.get(sizeY-1, x).set(Ex, Ey/c, B);
        }

        for(int y=1; y<sizeY-1; y++){
            Ex = -vPotentialField.get(y, 0).x;
            Ey = -0.5f*(PotentialField.get(y + 1, 0).x - PotentialField.get(y - 1, 0).x);

            Ex -= vPotentialField.get(y, 0).y;
            Ey -= vPotentialField.get(y, 0).z/c;

            B = vPotentialField.get(y, 0).z/c -
                    - 0.5f*(PotentialField.get(y + 1, 0).y - PotentialField.get(y - 1, 0).y);

            EMGravityField.get(y, 0).set(Ex/c, Ey, B);

            Ex = vPotentialField.get(y, 0).x;
            Ey = -0.5f*(PotentialField.get(y + 1, 0).x - PotentialField.get(y - 1, 0).x);

            Ex -= vPotentialField.get(y, 0).y;
            Ey -= vPotentialField.get(y, 0).z/c;

            B = -vPotentialField.get(y, 0).z/c -
                    - 0.5f*(PotentialField.get(y + 1, 0).y - PotentialField.get(y - 1, 0).y);

            EMGravityField.get(y, 0).set(Ex/c, Ey, B);
        }

        //Bottom-Left
        Ey = Ex = -vPotentialField.get(0, 0).x;

        Ex -= vPotentialField.get(0, 0).y;
        Ey -= vPotentialField.get(0, 0).z;

        B = vPotentialField.get(0,0).z
                - vPotentialField.get(0, 0).y;

        EMGravityField.get(0, 0).set(Ex, Ey, B).scl(1/c);

        //Bottom-Right
        Ex = vPotentialField.get(0, sizeX-1).x;
        Ey = -Ex;

        Ex -= vPotentialField.get(0, sizeX-1).y;
        Ey -= vPotentialField.get(0, sizeX-1).z;

        B = -vPotentialField.get(0,sizeX-1).z
                - vPotentialField.get(0, sizeX-1).y;

        EMGravityField.get(0, sizeX-1).set(Ex, Ey, B).scl(1/c);

        //Top-Left
        Ex = -vPotentialField.get(sizeY-1, 0).x;
        Ey = -Ex;

        Ex -= vPotentialField.get(sizeY-1, 0).y;
        Ey -= vPotentialField.get(sizeY-1, 0).z;

        B = vPotentialField.get(sizeY-1,0).z
                + vPotentialField.get(sizeY-1, 0).y;

        EMGravityField.get(sizeY-1, 0).set(Ex, Ey, B).scl(1/c);

        //Top-Right
        Ey = Ex = vPotentialField.get(sizeY-1, sizeX-1).x;

        Ex -= vPotentialField.get(sizeY-1, sizeX-1).y;
        Ey -= vPotentialField.get(sizeY-1, sizeX-1).z;

        B = -vPotentialField.get(sizeY-1, sizeX-1).z
                + vPotentialField.get(sizeY-1, sizeX-1).y;

        EMGravityField.get(sizeY-1, sizeX-1).set(Ex, Ey, B).scl(1/c);
    }

    public Vec2 Accel(int y, int x){
        Vec3 vF = EMGravityField.get(y, x);
        float Ex = vF.x,
              Ey = vF.y,
              B  = vF.z,
              xcomp = Tmp.v1.x, ycomp = Tmp.v1.y, xcomp2, ycomp2, dt = Time.delta*timescale;

        xcomp += (Ex + k*B*Tmp.v1.y)*dt;
        ycomp += (Ey - k*B*Tmp.v1.x)*dt;

        for(int i=0; i<15; i++){
            xcomp2 = Tmp.v1.x + (Ex + 0.5f*k*B*(Tmp.v1.y + ycomp)/c)*dt;
            ycomp2 = Tmp.v1.y + (Ey - 0.5f*k*B*(Tmp.v1.x + xcomp)/c)*dt;

            xcomp = xcomp2; ycomp = ycomp2;
        }

        return Tmp.v2.set(xcomp, ycomp).scl(dt);
    }

    public void update(){
        if(world.isGenerating()) return;
        Teams.TeamData D;
        Seq<Unit> U;
        Unit u;
        float x, y, m, w, h;
        int ix, iy;
        CurrentField.clear();
        for (Team t : Team.all) {
            D = t.data();
            if(D == null) continue;
            U = D.units;
            for(int i = 0; i < U.size; i++){
                u = U.get(i);
                x = u.x/tilesize;
                y = u.y/tilesize;
                ix = round(x); iy = round(y);
                if (ix < 1 || ix >= sizeX || iy < 1 || iy >= sizeY || Float.isNaN(x) || Float.isNaN(y)) continue;
                if(u.vel.len() >= maxSpeed) u.vel.setLength(maxSpeed);
                Tmp.v1.set(u.vel);

                m = G * u.mass()/ tilesize / tilesize / Mathf.sqrt(1f - Tmp.v1.len2()/csquared);
                w = x - ix; h = y - iy;

                CurrentField.get(iy-1, ix-1).add(csquared * m, Tmp.v1.x * m * (0.5f-w)*(0.5f-h), Tmp.v1.y * m * (0.5f-w)*(0.5f-h));
                CurrentField.get(iy-1, ix).add(csquared * m, Tmp.v1.x * m * (0.5f+w)*(0.5f-h), Tmp.v1.y * m * (0.5f+w)*(0.5f-h));
                CurrentField.get(iy, ix-1).add(csquared * m, Tmp.v1.x * m * (0.5f-w)*(0.5f+h), Tmp.v1.y * m * (0.5f-w)*(0.5f+h));
                CurrentField.get(iy, ix).add(csquared * m, Tmp.v1.x * m*(0.5f+w)*(0.5f+h), Tmp.v1.y * m * (0.5f+w)*(0.5f+h));
            }
        }

        float dt = Time.delta*timescale;

        //Yoshida Integrator

        BufferField.laplacian(PotentialField, vPotentialField).scl(stabilizer);
        PotentialField.addmul(BufferField.add(vPotentialField), c1*dt);
        BufferField.laplacian(PotentialField, vPotentialField);
        vPotentialField.addmul(BufferField.sub(CurrentField), d1*dt);

        BufferField.laplacian(PotentialField, vPotentialField).scl(stabilizer);
        PotentialField.addmul(BufferField.add(vPotentialField), c2*dt);
        BufferField.laplacian(PotentialField, vPotentialField);
        vPotentialField.addmul(BufferField.sub(CurrentField), d2*dt);

        BufferField.laplacian(PotentialField, vPotentialField).scl(stabilizer);
        PotentialField.addmul(BufferField.add(vPotentialField), c3*dt);
        BufferField.laplacian(PotentialField, vPotentialField);
        vPotentialField.addmul(BufferField.sub(CurrentField), d3*dt);

        BufferField.laplacian(PotentialField, vPotentialField).scl(stabilizer);
        PotentialField.addmul(BufferField.add(vPotentialField), c4*dt);

        //PotentialField.normalize();

        EMGravity();

        for (Team t : Team.all){
            D = t.data();
            if (D == null) continue;
            U = D.units;
            for (int i = 0; i < U.size; i++){
                u = U.get(i);
                x = u.x/tilesize;
                y = u.y/tilesize;
                ix = round(x); iy = round(y);
                w = x - ix; h = y - iy;

                if (ix < 1 || ix >= sizeX || iy < 1 || iy >= sizeY || Float.isNaN(x) || Float.isNaN(y)) continue;

                u.vel.add(Accel(iy - 1, ix - 1).scl((0.5f - w) * (0.5f - h)));
                u.vel.add(Accel(iy - 1, ix).scl((0.5f + w) * (0.5f - h)));
                u.vel.add(Accel(iy, ix - 1).scl((0.5f - w) * (0.5f + h)));
                u.vel.add(Accel(iy, ix).scl((0.5f + w) * (0.5f + h)));

                debug = debug.concat(u.vel.toString()).concat("\n");
            }
        }
    }

    @Override
    public int color(int i, int j){
        int color = 0xff, col = (int)(255f*Mathf.clamp(0.5f*EMGravityField.get(i, j).x/drawRadius + 0.5f));
        color |= col << 8;
        color |= col << 16;
        color |= col << 24;
        return color;
    }

    @Override
    public void write(DataOutput write) throws IOException{
        if (PotentialField == null || vPotentialField == null) return;

        write.writeInt(sizeX);
        write.writeInt(sizeY);

        PotentialField.eachTile(v -> {
            try{
                write.writeFloat(v.x);
                write.writeFloat(v.y);
                write.writeFloat(v.z);
            }catch(IOException e){
                throw new RuntimeException(e);
            }
        });
        vPotentialField.eachTile(v -> {
            try{
                write.writeFloat(v.x);
                write.writeFloat(v.y);
                write.writeFloat(v.z);
            }catch(IOException e){
                throw new RuntimeException(e);
            }
        });
    }

    @Override
    public void read(DataInput read) throws IOException{
        sizeX = read.readInt();
        sizeY = read.readInt();

        PotentialField = new GEMgrid(sizeY, sizeX);
        vPotentialField = new GEMgrid(sizeY, sizeX);
        BufferField = new GEMgrid(sizeY, sizeX);
        CurrentField = new GEMgrid(sizeY, sizeX);

        PotentialField.eachTile(v -> {
            try{
                v.x = read.readFloat();
                v.y = read.readFloat();
                v.z = read.readFloat();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
        vPotentialField.eachTile(v -> {
            try{
                v.x = read.readFloat();
                v.y = read.readFloat();
                v.z = read.readFloat();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
    }

    public static class GEMgrid{
        private final Vec3[] grid;
        private final Vec3 buffer = new Vec3();
        public int sizeY, sizeX;

        public GEMgrid(int sizeY, int sizeX) {
            grid = new Vec3[sizeX * sizeY];
            int totalSize = sizeX * sizeY;
            for (int i = 0; i < totalSize; i++){
                grid[i] = new Vec3();
            }
            this.sizeY = sizeY;
            this.sizeX = sizeX;
        }

        public final Vec3 get(int y, int x){
            return grid[x + y * sizeX];
        }

        public GEMgrid clear() {
            eachTile(v -> v.set(Vec3.Zero));
            return this;
        }

        /**
         * Note that the GEMgrid object that called this method will have its value be completely replaced.
         *
         * @param grid:  the grid that will be the values for the laplacian.
         * @param Bgrid: grid that will be used to impose the Neumann Boundary Condition: ∂(grid)/∂n = -(Bgrid)/c.
         */
        public GEMgrid laplacian(GEMgrid grid, GEMgrid Bgrid) {
            for (int i = 1; i < sizeY - 1; i++) {
                for (int j = 1; j < sizeX - 1; j++){
                    get(i, j).set(grid.get(i, j)).scl(-4f)
                            .add(grid.get(i, j - 1))
                            .add(grid.get(i, j + 1))
                            .add(grid.get(i - 1, j))
                            .add(grid.get(i + 1, j))
                            .scl(csquared);
                }
            }

            for (int i = 1; i < sizeY - 1; i++) { //X-edge
                buffer.set(Bgrid.get(i, 0)).scl(-2f / c);
                get(i, 0).set(grid.get(i, 0)).scl(-4f)
                        .add(grid.get(i, 1))
                        .add(grid.get(i, 1))
                        .add(grid.get(i - 1, 0))
                        .add(grid.get(i + 1, 0))
                        .add(buffer).scl(csquared);

                buffer.set(Bgrid.get(i, sizeX - 1)).scl(-2f / c);
                get(i, sizeX - 1).set(grid.get(i, sizeX - 1)).scl(-4f)
                        .add(grid.get(i, sizeX - 2))
                        .add(grid.get(i, sizeX - 2))
                        .add(grid.get(i - 1, sizeX - 1))
                        .add(grid.get(i + 1, sizeX - 1))
                        .add(buffer).scl(csquared);
            }


            for (int j = 1; j < sizeX - 1; j++) { //Y-edge
                buffer.set(Bgrid.get(0, j)).scl(-2f / c);
                get(0, j).set(grid.get(0, j)).scl(-4f)
                        .add(grid.get(1, j))
                        .add(grid.get(1, j))
                        .add(grid.get(0, j - 1))
                        .add(grid.get(0, j + 1))
                        .add(buffer).scl(csquared);

                buffer.set(Bgrid.get(sizeY - 1, j)).scl(-2f / c);
                get(sizeY - 1, j).set(grid.get(sizeY - 1, j)).scl(-4f)
                        .add(grid.get(sizeY - 2, j))
                        .add(grid.get(sizeY - 2, j))
                        .add(grid.get(sizeY - 1, j - 1))
                        .add(grid.get(sizeY - 1, j + 1))
                        .add(buffer).scl(csquared);
            }


            //Bottom-left
            get(0, 0).set(Bgrid.get(0, 0)).scl(-1f / c).sub(grid.get(0, 0)).scl(2f)
                    .add(grid.get(0, 1)).add(grid.get(1, 0)).scl(2f * csquared);
            //Bottom-right
            get(0, sizeX - 1).set(Bgrid.get(0, sizeX - 1)).scl(-1f / c).sub(grid.get(0, sizeX - 1)).scl(2f)
                    .add(grid.get(0, sizeX - 2)).add(grid.get(1, sizeX - 1)).scl(2f * csquared);
            //Top-left
            get(sizeY - 1, 0).set(Bgrid.get(sizeY - 1, 0)).scl(-1f / c).sub(grid.get(sizeY - 1, 0)).scl(2f)
                    .add(grid.get(sizeY - 1, 1)).add(grid.get(sizeY - 2, 0)).scl(2f * csquared);
            //Top-right
            get(sizeY - 1, sizeX - 1).set(Bgrid.get(sizeY - 1, sizeX - 1)).scl(-1f / c).sub(grid.get(sizeY - 1, sizeX - 1)).scl(2f)
                    .add(grid.get(sizeY - 1, sizeX - 2)).add(grid.get(sizeY - 2, sizeX - 1)).scl(2f * csquared);

            return this;
        }

        public GEMgrid set(GEMgrid grid) {
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).set(grid.get(i, j));
                }
            }
            return this;
        }

        public GEMgrid add(GEMgrid grid) {
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).add(grid.get(i, j));
                }
            }
            return this;
        }

        public GEMgrid sub(GEMgrid grid) {
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).sub(grid.get(i, j));
                }
            }
            return this;
        }

        public GEMgrid sub(Vec3 v){
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).sub(v);
                }
            }
            return this;
        }

        public GEMgrid normalize(){
            buffer.set(Vec3.Zero);
            final float scale = 1f/grid.length;
            eachTile(v -> buffer.add(v, scale));
            sub(buffer);
            return this;
        }

        public GEMgrid scl(float s) {
            eachTile(v -> v.scl(s));
            return this;
        }

        public GEMgrid addmul(GEMgrid grid, float s){
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).add(grid.get(i, j), s);
                }
            }
            return this;
        }

        public GEMgrid eachTile(Cons<Vec3> cons) {
            for (Vec3 v : grid) {
                cons.get(v);
            }
            return this;
        }
    }
}
